#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <iomanip>

using namespace std;

// Parameter wie zuvor
const double v = 0.2;
const double L = 1.0;
const double dt = 0.01;
const double x0 = 0.5;
const double initialPulseArea = 1.0;
const double D3 = 0.0005;  // Diffusionskoeffizient für Stoff 3 (kritischster Fall)

// Analytische Lösung (für Stoff 3 mit D3)
double analytischeLoesung(double D, double t, double x) {
    if (t == 0) {
        return (fabs(x - x0) < 1e-12 ? initialPulseArea/(/*dx*/1.0) : 0.0); 
        // Hinweis: dx wird hier nicht verwendet, da wir bei t=0 nur den Dirac-Puls qualitativ darstellen.
    }
    double denom = sqrt(4.0 * M_PI * D * t);
    double exponent = - pow(x - x0 - v*t, 2.0) / (4.0 * D * t);
    return 1.0/denom * exp(exponent);
}

int main() {
    // Verschiedene räumliche Schrittweiten zur Untersuchung
    double dxValues[] = {0.02, 0.01, 0.005};
    int nCases = 3;
    // Ausgabe der numerischen Peclet-Zahlen für jede Diskretisierung
    cout << "Numerische Peclet-Zahlen für verschiedene dx (Stoff 3):\n";
    for (int idx = 0; idx < nCases; ++idx) {
        double dx = dxValues[idx];
        double Pe = v * dx / (2.0 * D3);
        cout << "dx = " << dx << ": Pe_num = " << Pe << "\n";
    }
    cout << endl;
    // Simulation für jede dx-Stufe
    for (int idx = 0; idx < nCases; ++idx) {
        double dx = dxValues[idx];
        int N = int(L/dx) + 1;
        int steps = int(1.0/dt);
        // Vektoren für numerische Lösung (für Zentral und Upwind getrennt)
        vector<double> uCentral(N), uCentral_new(N);
        vector<double> uUpwind(N), uUpwind_new(N);
        // Anfangsbedingungen: Dirac-Puls an x0 = 0.5
        int i0 = (int) round(x0 / dx);
        for (int i = 0; i < N; ++i) {
            if (i == i0) {
                uCentral[i] = initialPulseArea / dx;
                uUpwind[i] = initialPulseArea / dx;
            } else {
                uCentral[i] = 0.0;
                uUpwind[i] = 0.0;
            }
        }
        // Zeitschleife bis t=1.0
        double t = 0.0;
        for (int step = 1; step <= steps; ++step) {
            // Randbedingungen: Dirichlet an x=0 und x=1 (beide 0 für unseren Fall)
            uCentral[0] = 0.0;
            uCentral[N-1] = 0.0;
            uUpwind[0] = 0.0;
            uUpwind[N-1] = 0.0;
            // RHS für impliziten Schritt berechnen:
            // Zentraldifferenzen (Konvektion):
            vector<double> rhsCentral(N);
            // Upwind (Konvektion):
            vector<double> rhsUpwind(N);
            rhsCentral[0] = 0.0;
            rhsCentral[N-1] = 0.0;
            rhsUpwind[0] = 0.0;
            rhsUpwind[N-1] = 0.0;
            for (int i = 1; i < N-1; ++i) {
                // Zentraldifferenzen Konvektionsterm:
                double convCentral = v * (uCentral[i+1] - uCentral[i-1]) / (2.0 * dx);
                rhsCentral[i] = uCentral[i] - dt * convCentral;
                // Upwind-Schema Konvektionsterm (v>0: Rueckwärtsdifferenz):
                double convUp = v * (uUpwind[i] - uUpwind[i-1]) / dx;
                rhsUpwind[i] = uUpwind[i] - dt * convUp;
            }
            // Impliziter Diffusionsschritt (D = D3 konstant):
            double r = D3 * dt / (dx * dx);
            int n_interior = N - 2;
            // Tridiagonale Koeffizienten (gleich für beide Schemata, da Diffusionsteil identisch)
            vector<double> a(n_interior + 1), b(n_interior + 1), c(n_interior + 1);
            for (int innerIndex = 1; innerIndex <= n_interior; ++innerIndex) {
                if (innerIndex == 1) a[innerIndex] = 0.0;
                else a[innerIndex] = -r;
                if (innerIndex == n_interior) c[innerIndex] = 0.0;
                else c[innerIndex] = -r;
                b[innerIndex] = 1.0 + 2.0 * r;
            }
            // Thomas-Algorithmus für Zentraldif.-Lösung
            vector<double> d_c(n_interior + 1);
            for (int innerIndex = 1; innerIndex <= n_interior; ++innerIndex) {
                d_c[innerIndex] = rhsCentral[innerIndex];
            }
            for (int k = 2; k <= n_interior; ++k) {
                double m = a[k] / b[k-1];
                b[k] -= m * c[k-1];
                d_c[k] -= m * d_c[k-1];
            }
            vector<double> x_c(n_interior + 1);
            x_c[n_interior] = d_c[n_interior] / b[n_interior];
            for (int k = n_interior - 1; k >= 1; --k) {
                x_c[k] = (d_c[k] - c[k] * x_c[k+1]) / b[k];
            }
            // Thomas-Algorithmus für Upwind-Lösung (Koeff. b und a,c identisch, aber d unterschiedlich)
            vector<double> d_u(n_interior + 1);
            for (int innerIndex = 1; innerIndex <= n_interior; ++innerIndex) {
                d_u[innerIndex] = rhsUpwind[innerIndex];
            }
            for (int k = 2; k <= n_interior; ++k) {
                double m = a[k] / (/* Achtung: b wurde oben modifiziert, daher neu initialisieren für Upwind */ (1.0 + 2.0 * r) - (k>2 ? (a[k] / b[k-1]) * c[k-1] : 0.0));
                // Hinweis: Um die Elimination getrennt zu behandeln, führen wir sie erneut durch:
            }
            // **Vereinfachung**: Führe die Elimination für Upwind separat noch einmal vollständig durch
        }
    }
    return 0;
}
