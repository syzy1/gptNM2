#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>  // für ConjugateGradient
#include "Library.h"                     // alle generate*, solve*, writeVectorAsMatrix, gnuplot_splot

using namespace std;
using namespace Eigen;
using Clock = chrono::high_resolution_clock;

int main(){
    vector<int> ns = {10, 25, 40, 60};
    double fVal = 1.0, gVal = 0.0;
    int m = ns.size();

    // Report-Datei
    ofstream out("solver_comparison.txt");
    out << "n  timeLU[s]  timeCG[s]  itersCG\n";

    for(int n : ns){
        // 1) System aufbauen
        SparseMatrix<double> A = generatePoissonMatrix2D(n);
        VectorXd b             = generateRHS2D(n, fVal, gVal);

        // 2) Direkter Solver (SparseLU) messen
        auto t0 = Clock::now();
        VectorXd xLU = solveSparseLU(A, b);
        auto t1 = Clock::now();
        double timeLU = chrono::duration<double>(t1 - t0).count();

        // Lösung speichern
        writeVectorAsMatrix("solution_LU_" + to_string(n) + ".txt",
                            xLU, n-1, n-1);

        // 3) Iterativer Solver (CG) messen und Iterationen abfragen
        ConjugateGradient<SparseMatrix<double>, Lower|Upper> cg;
        auto t2 = Clock::now();
        cg.compute(A);
        VectorXd xCG = cg.solve(b);
        auto t3 = Clock::now();
        double timeCG = chrono::duration<double>(t3 - t2).count();
        int iterCG = cg.iterations();

        // CG-Lösung speichern
        writeVectorAsMatrix("solution_CG_" + to_string(n) + ".txt",
                            xCG, n-1, n-1);

        // 4) In Report schreiben
        out << n << "  "
            << timeLU << "  "
            << timeCG << "  "
            << iterCG << "\n";
    }

    out.close();
    cout << "Fertig. Daten in solver_comparison.txt\n";

    // 5) Plots erstellen
    //gnuplot_splot("solver_comparison.txt");
    return 0;
}
