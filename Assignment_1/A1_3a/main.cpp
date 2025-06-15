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

    // --- Anfangswerte
    VectorXd y0(3);
    y0 << 1.0, 0.0, 0.0;

    // --- Solver mit Richardson-Extrapolation
    double tol = 1e-6;
    int kmax = 1000;
    Euler solver(h, tol, kmax, ReactionODE);

    // semiimplizit_do führt interne Halbierung und Extrapolation durch
    MatrixXd Y = solver.semiimplizit_do(y0, tVec);
    string outFile = "SemiImplizit_Richardson.txt";
    writeVectorMatrix(tVec, Y, outFile);
    cout << "Semi-Implizit (Richardson) in '"<< outFile <<"' gespeichert.\n";

    // --- Gnuplot
    ofstream gp("plot_semiimplizit.gp");
    gp << "reset\n"
       << "set title 'Semi-Implizites Euler + Richardson'\n"
       << "set xlabel 't [s]'\n"
       << "set ylabel 'c [mol/L]'\n"
       << "plot '"<< outFile <<"' using 1:2 with lines title 'cA', "
          "'" << outFile <<"' using 1:3 with lines title 'cB', "
          "'" << outFile <<"' using 1:4 with lines title 'cC'\n";
    gp.close();
    system("gnuplot plot_semiimplizit.gp -");

    return 0;
}
