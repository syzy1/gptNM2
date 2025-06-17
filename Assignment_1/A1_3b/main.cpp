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
static double k1 = 1.0;
static double k2 = 10.0;

// Definition des ODE-Systems
void ReactionODE(VectorXd &y, VectorXd &F) {
    double cA = y(0), cB = y(1), cC = y(2);
    F(0) = -k1 * cA;
    F(1) =  k1 * cA - k2 * cB;
    F(2) =  k2 * cB;
}

int main() {
    // 1) Zeitgitter
    double h    = 0.05, t0 = 0.0, tEnd = 20.0;
    int    n    = int((tEnd - t0)/h) + 1;
    VectorXd t  = VectorXd::LinSpaced(n, t0, tEnd);
    int    n2   = 2*(n-1) + 1;
    VectorXd t2 = VectorXd::LinSpaced(n2, t0, tEnd);

    // 2) Anfangswerte
    VectorXd y0(3);
    y0 << 1.0, 0.0, 0.0;

    // 3) Solver initialisieren
    double tol  = 1e-6;
    int    kmax = 1000;
    Euler solverLow  (h,     tol, kmax, ReactionODE);
    Euler solverHigh (h/2.0, tol, kmax, ReactionODE);

    // 4) Lösungen berechnen
    MatrixXd Y_low       = solverLow .semiimplizit(y0, t);
    MatrixXd Y_high_full = solverHigh.semiimplizit(y0, t2);
    MatrixXd Y_high(n, y0.size());
    for (int i = 0; i < n; ++i)
        Y_high.row(i) = Y_high_full.row(2*i);
    MatrixXd Y_extr = solverLow.semiimplizit_do(y0, t);

    // 5) Fehlerschätzungen
    VectorXd E_low(n), E_high(n);
    for (int i = 0; i < n; ++i) {
        E_low(i)  = (Y_extr.row(i) - Y_low.row(i)).norm();
        E_high(i) = (Y_extr.row(i) - Y_high.row(i)).norm();
    }

    // 6) Ergebnisse in Datei
    ofstream fout("Ergebnisse.txt");
    fout << "# t[s]\t"
         << "cA_low\tcB_low\tcC_low\t"
         << "cA_high\tcB_high\tcC_high\t"
         << "cA_extr\tcB_extr\tcC_extr\t"
         << "Err_low\tErr_high\n";
    for (int i = 0; i < n; ++i) {
        fout << t(i)           << "\t"
             << Y_low(i,0)     << "\t" << Y_low(i,1)     << "\t" << Y_low(i,2)     << "\t"
             << Y_high(i,0)    << "\t" << Y_high(i,1)    << "\t" << Y_high(i,2)    << "\t"
             << Y_extr(i,0)    << "\t" << Y_extr(i,1)    << "\t" << Y_extr(i,2)    << "\t"
             << E_low(i)       << "\t" << E_high(i)      << "\n";
    }
    fout.close();

    // 7a) Plot cA bis t=3s
    ofstream gpA("plot_cA_zoom.gp");
    gpA << "reset\n"
        << "set title 'Zeitverlauf von cA (bis 3s)'\n"
        << "set xlabel 't [s]'\n"
        << "set ylabel 'cA [mol/L]'\n"
        << "set xrange [0:3]\n"
        << "plot 'Ergebnisse.txt' using 1:2 with lines title 'cAlow', \\\n"
        << "     '' using 1:5 with lines title 'cAhigh', \\\n"
        << "     '' using 1:8 with lines title 'cAextr'\n";
    gpA.close();
    system("gnuplot -persist plot_cA_zoom.gp");

    // 7b) Plot cB bis t=3s
    ofstream gpB("plot_cB_zoom.gp");
    gpB << "reset\n"
        << "set title 'Zeitverlauf von cB (bis 3s)'\n"
        << "set xlabel 't [s]'\n"
        << "set ylabel 'cB [mol/L]'\n"
        << "set xrange [0:3]\n"
        << "plot 'Ergebnisse.txt' using 1:3 with lines title 'cBlow', \\\n"
        << "     '' using 1:6 with lines title 'cBhigh', \\\n"
        << "     '' using 1:9 with lines title 'cBextr'\n";
    gpB.close();
    system("gnuplot -persist plot_cB_zoom.gp");

    // 7c) Plot cC bis t=3s
    ofstream gpC("plot_cC_zoom.gp");
    gpC << "reset\n"
        << "set title 'Zeitverlauf von cC (bis 3s)'\n"
        << "set xlabel 't [s]'\n"
        << "set ylabel 'cC [mol/L]'\n"
        << "set xrange [0:3]\n"
        << "plot 'Ergebnisse.txt' using 1:4 with lines title 'cClow', \\\n"
        << "     '' using 1:7 with lines title 'cChigh', \\\n"
        << "     '' using 1:10 with lines title 'cCextr'\n";
    gpC.close();
    system("gnuplot -persist plot_cC_zoom.gp");

    // 7d) Fehlerplot (unverändert)
    ofstream gpE("plot_error.gp");
    gpE << "reset\n"
        << "set title 'Fehler gegen low/extr und high/extr'\n"
        << "set xlabel 't [s]'\n"
        << "set ylabel 'Error'\n"
        << "plot 'Ergebnisse.txt' using 1:11 with lines title 'Errlow', \\\n"
        << "     '' using 1:12 with lines title 'Errhigh'\n";
    gpE.close();
    system("gnuplot -persist plot_error.gp");

    return 0;
}
