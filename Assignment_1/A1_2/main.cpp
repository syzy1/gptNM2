#include <iostream>
#include <fstream>
#include <string>
#include "Library.h"
#include "Euler.h"
#include "Methods.h"

using namespace std;
using namespace Eigen;

static double k1 = 1.0;
static double k2 = 10.0;

void ReactionODE(VectorXd &y, VectorXd &F) {
    double cA = y(0), cB = y(1), cC = y(2);
    F(0) = -k1*cA;
    F(1) =  k1*cA - k2*cB;
    F(2) =  k2*cB;
}

int main() {
    // --- Zeitgitter
    double h = 0.05, t0 = 0.0, tEnd = 20.0;
    int nSteps = int((tEnd - t0)/h) + 1;
    VectorXd tVec(nSteps);
    for(int i=0; i<nSteps; ++i) tVec(i) = t0 + i*h;

    VectorXd y0(3);
    y0 << 1.0, 0.0, 0.0;

    double tol = 1e-6;
    int kmax = 1000;
    Euler solver(h, tol, kmax, ReactionODE);

    // --- Auswahlmenü
    cout << "Wählen Sie das Verfahren:\n"
         << " 1) Explizites Euler\n"
         << " 2) Implizites Euler\n"
         << "Eingabe: ";
    int wahl; 
    cin >> wahl;

    MatrixXd Y;
    string methodName, outFile;
    if (wahl == 1) {
        Y = solver.explizit(y0, tVec);
        methodName = "Explizit";
        outFile    = "Wahl_Explizit.txt";
    } else {
        Y = solver.implizit(y0, tVec);
        methodName = "Implizit";
        outFile    = "Wahl_Implizit.txt";
    }

    writeVectorMatrix(tVec, Y, outFile);
    cout << "Ergebnisse ("<< methodName <<") in '"<< outFile <<"' gespeichert.\n";

    // --- Gnuplot
    ofstream gp("plot_wahl.gp");
    gp << "reset\n"
       << "set title 'Euler-Verfahren ("<< methodName <<")'\n"
       << "set xlabel 't [s]'\n"
       << "set ylabel 'c [mol/L]'\n"
       << "plot '"<< outFile <<"' using 1:2 with lines title 'cA', "
          "'" << outFile <<"' using 1:3 with lines title 'cB', "
          "'" << outFile <<"' using 1:4 with lines title 'cC'\n";
    gp.close();
    system("gnuplot plot_wahl.gp -");

    return 0;
}
