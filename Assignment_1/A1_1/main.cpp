#include <iostream>      // Ein-/Ausgabe
#include <fstream>       // Datei-Operationen
#include <string>
#include "Library.h"     // liefert VectorXd, MatrixXd, etc.
#include "Euler.h"       // explizites Euler-Verfahren
#include "Methods.h"     // writeVectorMatrix

using namespace std;
using namespace Eigen;

// Reaktionsgeschwindigkeit
static double k1 = 1.0;
static double k2 = 10.0;

// rechte Seite der DGL y' = f(t,y)
void ReactionODE(VectorXd &y, VectorXd &F) {
    double cA = y(0), cB = y(1), cC = y(2);
    F(0) = -k1*cA;
    F(1) =  k1*cA - k2*cB;
    F(2) =  k2*cB;
}

int main() {
    // --- Zeitgitter
    double h = 0.05;
    double t0 = 0.0, tEnd = 20.0;
    int nSteps = int((tEnd - t0)/h) + 1;
    VectorXd tVec(nSteps);
    for(int i=0; i<nSteps; ++i)
        tVec(i) = t0 + i*h;

    // --- Anfangswerte [cA, cB, cC]
    VectorXd y0(3);
    y0 << 1.0, 0.0, 0.0;

    // --- Solver einrichten
    double tol = 1e-6;
    int kmax = 1000;
    Euler solver(h, tol, kmax, ReactionODE);

    // --- Rechnung: explizites Euler
    MatrixXd Y = solver.explizit(y0, tVec);
    string outFile = "Explizit_Ergebnisse.txt";
    writeVectorMatrix(tVec, Y, outFile);
    cout << "Ergebnisse in '"<< outFile <<"' gespeichert.\n";

    // --- Gnuplot-Skript erzeugen und ausführen
    ofstream gp("plot_explizit.gp");
    gp
      << "reset\n"
      << "set title 'Explizites Euler-Verfahren'\n"
      << "set xlabel 't [s]'\n"
      << "set ylabel 'c [mol/L]'\n"
      << "plot '"<< outFile <<"' using 1:2 with lines title 'cA', "
         "'" << outFile <<"' using 1:3 with lines title 'cB', "
         "'" << outFile <<"' using 1:4 with lines title 'cC'\n";
    gp.close();
    system("gnuplot plot_explizit.gp -");

    return 0;
}
