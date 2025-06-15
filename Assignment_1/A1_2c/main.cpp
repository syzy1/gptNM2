#include <iostream>    // Für Ein- und Ausgabe
#include <fstream>     // Für Datei-Operationen
#include <cmath>
#include <string>
#include "Library.h"   // Enthält ggf. Newton, ReactionODE, etc.
#include "Euler.h"
#include "Methods.h"

using namespace std;
using namespace Eigen;

// Globale Reaktionskonstanten (werden in ReactionODE genutzt)
static double k1 = 1.0;
static double k2 = 1000000;  // k2 bleibt 10^6

// Reaktions-ODE: y=[cA, cB, cC], F= [dcA/dt, dcB/dt, dcC/dt].
void ReactionODE(VectorXd &y, VectorXd &F)
{
    double cA = y(0);
    double cB = y(1);
    double cC = y(2);

    F(0) = -k1 * cA;
    F(1) =  k1 * cA - k2 * cB;
    F(2) =  k2 * cB;
}

void runEulerMethods(double h, const string &h_label) {
    double t_start = 0.0;
    double t_end   = 1e-4;  // Bereich 0 ≤ t ≤ 10^-4 s
    int    nSteps  = static_cast<int>((t_end - t_start)/h) + 1;

    VectorXd tVec(nSteps);
    for(int i = 0; i < nSteps; i++) {
        tVec(i) = t_start + i * h;
    }

    // Anfangswerte: cA=1, cB=0, cC=0
    VectorXd y0(3);
    y0 << 1.0, 0.0, 0.0;

    // Euler-Objekt
    double tol = 1e-6;
    int    kmax = 1000;
    Euler solver(h, tol, kmax, ReactionODE);

    // Explizites Euler
    MatrixXd solution_explicit = solver.explizit(y0, tVec);
    string   outfile_explizit  = "Ergebnisse_Explizit_h=" + h_label + ".txt";
    writeVectorMatrix(tVec, solution_explicit, outfile_explizit);
    cout << "Explizites Euler mit h=" << h << " gespeichert in " << outfile_explizit << endl;

    // Implizites Euler
    MatrixXd solution_implicit = solver.implizit(y0, tVec);
    string   outfile_implicit  = "Ergebnisse_Implizit_h=" + h_label + ".txt";
    writeVectorMatrix(tVec, solution_implicit, outfile_implicit);
    cout << "Implizites Euler mit h=" << h << " gespeichert in " << outfile_implicit << endl;
}

void generateGnuplotScript() {
    ofstream gp("plot_commands.gp");
    gp << "reset\n";
    gp << "set terminal qt\n";  // Fenster-Plot
    gp << "set xlabel 't [s]'\n";
    gp << "set ylabel 'Konzentration [mol/L]'\n";
    gp << "set logscale y\n";
    gp << "set yrange [1e-8:1.2]\n";

    // Plot für alle 6 Ergebnisse
    vector<string> methods   = {"Explizit", "Implizit"};
    vector<string> stepSizes = {"2e-6", "1e-6", "5e-7"};

    for (const auto& method : methods) {
        for (const auto& h_label : stepSizes) {
            string filename = "Ergebnisse_" + method + "_h=" + h_label + ".txt";

            gp << "set title '" << method << " Euler, h=" << h_label << "'\n";
            gp << "plot '" << filename << "' using 1:2 with lines title 'cA', ";
            gp << "'" << filename << "' using 1:3 with lines title 'cB', ";
            gp << "'" << filename << "' using 1:4 with lines title 'cC'\n";
            gp << "pause mouse close\n"; // Stoppe nach jedem Plot
        }
    }

    gp.close();
}

int main() {
    cout << "Starte Berechnungen für verschiedene Schrittweiten...\n";

    runEulerMethods(2e-6, "2e-6");
    runEulerMethods(1e-6, "1e-6");
    runEulerMethods(5e-7, "5e-7");

    cout << "Alle Ergebnisse wurden gespeichert.\n";

    // Gnuplot-Skripte generieren
    generateGnuplotScript();

    // Mit -persist, damit alle Fenster nacheinander korrekt dargestellt werden
    system("gnuplot -persist plot_commands.gp");

    return 0;
}
