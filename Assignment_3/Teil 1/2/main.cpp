#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <iomanip>
#include "Library.h"    // Professor’s Bibliothek

using namespace std;

// physikalische Parameter
const double v   = 0.2;    // Strömungsgeschwindigkeit
const double L   = 1.0;    // Länge des Gebietes
const double dx  = 0.02;   // Zellbreite Δx
const double dt  = 0.01;   // Zeitschritt Δt
const int    N   = int(L/dx) + 1;  // Anzahl Gitterpunkte
const double x0  = 0.5;    // Puls-Mittelpunkt
const double A0  = 1.0;    // Fläche des Anfangspulses

// analytische Lösung für Vergleich
double analytischeLoesung(double D, double t, double x) {
    if (t == 0.0) {
        // bei t=0: Dirac-Puls an x0
        return (fabs(x - x0) < 1e-12 ? A0/dx : 0.0);
    }
    const double pi = acos(-1.0);                      // Kreiszahl π
    double denom   = sqrt(4.0 * pi * D * t);           // Normalisierung
    double arg     = -pow(x - x0 - v*t, 2.0) 
                     / (4.0 * D * t);                  // Exponenten-Term
    return (1.0/denom) * exp(arg);                     // Gauß-Funktion
}

int main() {
    int i0 = int(round(x0/dx));    // Index der Zelle mit Puls

    // Lösungsvektoren und Fluss-Array
    vector<double> uOld(N),       // alte Zeitschicht
                   uNew(N),       // neue Zeitschicht
                   flux(N+1);     // Fluss an Zellgrenzen

    // Schleife über drei Stoffe
    for (int j = 1; j <= 3; ++j) {
        double D = (j==1 ? 0.01      // Diffusionskoeffizient Stoff 1
                  : (j==2 ? 0.005    // Diffusionskoeffizient Stoff 2
                           : 0.0005)); // Diffusionskoeffizient Stoff 3
        double Pe = v * dx / D;       // Péclet-Zahl
        cout << "Stoff " << j << "  Pe = " << Pe << "\n";

        // zwei Schemata: Linear (zentral) und Upwind
        for (bool useUpwind : {false, true}) {
            string scheme = useUpwind 
                            ? "Upwind" 
                            : "Linear";             // Schema-Name

            // Ausgabedateien für t=0.5 und t=1.0
            ofstream f05("Stoff" + to_string(j) 
                         + "_" + scheme + "_t05.txt");
            ofstream f10("Stoff" + to_string(j) 
                         + "_" + scheme + "_t10.txt");
            f05 << scientific << setprecision(6);
            f10 << scientific << setprecision(6);

            // Anfangsbedingung: Dirac-Puls in Zelle i0
            for (int i = 0; i < N; ++i)
                uOld[i] = (i==i0 ? A0/dx : 0.0);

            int totalSteps = int(1.0/dt);  // Anzahl Zeitschritte
            for (int step = 1; step <= totalSteps; ++step) {
                double t = step * dt;       // aktuelle Zeit

                // Dirichlet-Rand: u=0 an Domain-Enden
                uOld[0] = uOld[N-1] = 0.0;

                // Fluss φ_{i+½} = v·u_face − D·(u_R − u_L)/Δx
                flux[0] = 0.0;              // linker Rand
                flux[N] = 0.0;              // rechter Rand
                for (int f = 1; f < N; ++f) {
                    double uL = uOld[f-1];  // linker Zellwert
                    double uR = uOld[f];    // rechter Zellwert
                    // Wahl des Mittelwerts für Konvektion
                    double uFace = useUpwind
                                   ? (v >= 0 ? uL : uR) 
                                   : 0.5*(uL + uR);
                    double diffFlux = -D * (uR - uL) / dx; // Diffusion
                    flux[f] = v * uFace + diffFlux;        // Gesamtfluss
                }

                // Finite-Volumen-Update
                for (int i = 1; i < N-1; ++i) {
                    uNew[i] = uOld[i]
                              - (dt/dx) 
                              * (flux[i+1] - flux[i]); // Massenbilanz
                }
                uNew[0] = uNew[N-1] = 0.0; // Ränder erneut setzen

                // Ausgabe bei t=0.5 und t=1.0
                if (step == int(0.5/dt)) {
                    for (int i = 0; i < N; ++i) {
                        double x = i * dx;    // Ortskoordinate
                        f05 << x << "\t"
                             << uNew[i] << "\t"
                             << analytischeLoesung(D, t, x)
                             << "\n";
                    }
                }
                if (step == totalSteps) {
                    for (int i = 0; i < N; ++i) {
                        double x = i * dx;
                        f10 << x << "\t"
                             << uNew[i] << "\t"
                             << analytischeLoesung(D, t, x)
                             << "\n";
                    }
                }

                uOld.swap(uNew); // für nächsten Schritt tauschen
            }

            f05.close();
            f10.close();
            cout << "  Schema " << scheme 
                 << " abgeschlossen.\n";
        }
    }

    return 0;
}
