#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "Library.h"   // liefert VectorXd, MatrixXd, solve‐Funktionen
#include "Euler.h"
#include "Methods.h"   // für jac_cd, wenn du direkt die Jacobi‐Matrix willst

using namespace std;
using namespace Eigen;

// globale Reaktions‐Konstanten, werden in ReactionODE genutzt
static double k1 = 1.0;
static double k2 = 10.0;

// RHS der ODE y' = f(y)
void ReactionODE(VectorXd &y, VectorXd &F) {
    double cA = y(0), cB = y(1), cC = y(2);
    F(0) = -k1 * cA;
    F(1) =  k1 * cA - k2 * cB;
    F(2) =       k2 * cB;
}

// Hilfsroutine, Zahlen “schön” als string
string formatZahl(double x) {
    if (fabs(x - 1e6) < 1e-9)   return "1e6";
    if (fabs(x + 1e6) < 1e-9)   return "-1e6";
    double r = round(x);
    if (fabs(x - r) < 1e-12)    return to_string((long)r);
    ostringstream oss; oss << x;
    return oss.str();
}

int main() {
    // Anfangszustand y0 = [1,0,0]
    VectorXd y0(3);
    y0 << 1.0, 0.0, 0.0;

    // Ausgabedatei
    ofstream out("Steifigkeit_Ergebnisse.txt");
    if (!out) {
        cerr << "Fehler: konnte Datei nicht öffnen\n";
        return 1;
    }
    out << "=== Steifigkeits-Analyse ===\n\n";

    // Lambda für jeden Fall
    auto doCase = [&](double K1, double K2, int id) {
        // Konstanten setzen
        k1 = K1;  k2 = K2;

        // Euler-Objekt für Jacobi-Differenzen
        double hJac   = 0.05;    // Schrittweite für jac_cd (lineare ODE → exakt)
        double tol    = 1e-6;    // wird hier nicht genutzt
        int    maxIt  = 10;      // auch nicht genutzt
        Euler solver(hJac, tol, maxIt, ReactionODE);

        // Jacobi-Matrix & Eigenwerte (reell)
        MatrixXd J = jac_cd(ReactionODE, y0, hJac);
        EigenSolver<MatrixXd> es(J);
        VectorXd evalsRe = es.eigenvalues().real();

        // Steifigkeitsquotient
        double SR = solver.Steifigkeit(y0);

        // Klassifikation
        string cls = (SR < 500.0)
                     ? "keine steife DGL"
                     : (SR < 1e5
                        ? "steife DGL"
                        : "sehr steife DGL");

        // Ausgabe
        out << "Fall " << id << ": k1=" << formatZahl(K1)
            << ", k2=" << formatZahl(K2) << "\n";
        out << "Eigenwerte:  ";
        for (int i = 0; i < evalsRe.size(); ++i) {
            out << formatZahl(evalsRe[i])
                << (i + 1 < evalsRe.size() ? "   " : "\n");
        }

        // max/min (absolut) für SR‐Formel
        double vmax = *max_element(evalsRe.data(),
                                   evalsRe.data()+evalsRe.size(),
                                   [](double a,double b){return fabs(a)<fabs(b);});
        // vmin: kleinstes |eig| > 0
        vector<double> nonZero;
        for (int i=0; i<evalsRe.size(); ++i)
            if (fabs(evalsRe[i])>1e-12) nonZero.push_back(fabs(evalsRe[i]));
        double vmin = *min_element(nonZero.begin(), nonZero.end());

        out << "SR = |" << formatZahl(vmax)
            << "|/|"  << formatZahl(vmin)
            << "| = " << formatZahl(SR) << "\n";
        out << "→ " << cls << "\n\n";
    };

    // beide Fälle
    doCase(1.0,   10.0, 1);
    doCase(1.0, 1e6,    2);

    out.close();
    cout << "Fertig: Steifigkeit_Ergebnisse.txt erstellt\n";
    return 0;
}
