#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "Library.h"
#include "Euler.h"
#include "Methods.h"

using namespace std;
using namespace Eigen;

// Reaktionsparameter
static double k1 = 1.0, k2 = 10.0;

// ODE-System
void ReactionODE(VectorXd &y, VectorXd &F) {
    double cA=y(0), cB=y(1), cC=y(2);
    F(0) = -k1*cA;
    F(1) =  k1*cA - k2*cB;
    F(2) =  k2*cB;
}

int main() {
    // 1) Zeitschritte
    double h    = 0.05;
    double t0   = 0.0, tEnd = 20.0;
    int    n    = int((tEnd - t0)/h) + 1;
    VectorXd t  = VectorXd::LinSpaced(n, t0, tEnd);

    // Anfangswerte
    VectorXd y0(3); 
    y0 << 1.0, 0.0, 0.0;

    // Solver
    double tol = 1e-6;
    int    kmax= 1000;
    Euler solverLow (h,     tol, kmax, ReactionODE);
    Euler solverHigh(h/2.0, tol, kmax, ReactionODE);

    // 2) Lösungen berechnen
    MatrixXd Y_low  = solverLow.semiimplizit(y0, t);
    int      n2     = 2*(n-1) + 1;
    VectorXd t2     = VectorXd::LinSpaced(n2, t0, tEnd);
    MatrixXd Y_high = solverHigh.semiimplizit(y0, t2);

    // 3) Fehlerschätzung und Kurz-Check auf Konsole
    vector<double> err(n);
    cout << "Erste 5 Fehlerwerte (t_i , E_i):\n";
    for(int i=0; i<n; ++i) {
        err[i] = (Y_high.row(2*i) - Y_low.row(i)).norm();
        if(i<5)
            cout << "  " << t(i) << "  ,  " << err[i] << "\n";
    }

    // 4) Datei-Ausgabe
    ofstream out("Fehler_Ergebnis.txt");
    out << "# t  cA_high  cB_high  cC_high  err_norm\n";
    for(int i=0; i<n; ++i) {
        out << t(i) << " "
            << Y_high(2*i,0) << " "
            << Y_high(2*i,1) << " "
            << Y_high(2*i,2) << " "
            << err[i] << "\n";
    }
    out.close();
    cout << "Vollständige Fehleranalyse in 'Fehler_Ergebnis.txt'.\n";

    // 5) Gnuplot-Skript für zwei Fenster erzeugen
    ofstream gp("plot_fehler.gp");
    gp 
      << "set term qt 0 title 'Konzentrationen (semi-implizit)'\n"
      << "plot 'Fehler_Ergebnis.txt' using 1:2 with lines title 'cA', \\\n"
      << "     '' using 1:3 with lines title 'cB', \\\n"
      << "     '' using 1:4 with lines title 'cC'\n\n"
      << "set term qt 1 title 'Fehler E_i = ||Y_high−Y_low||'\n"
      << "plot 'Fehler_Ergebnis.txt' using 1:5 with lines title 'Fehler'\n";
    gp.close();

    // 6) Gnuplot aufrufen und beide Fenster offen halten
    system("gnuplot -persist plot_fehler.gp");

    return 0;
}
