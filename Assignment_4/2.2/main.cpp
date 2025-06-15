#include <iostream>
#include <fstream>
#include "Methods.h"
#include "Library.h"  // enthält generatePoissonMatrix2D, generateRHS2D, solveSparseLU, solveSparseCG, writeVectorAsMatrix, gnuplot_splot
using namespace std;
using namespace Eigen;

int main(){
    // 1) n, f, g einlesen (Datei oder Konsole)
    int n; double f, g;
    ifstream fin("poisson_data.txt");
    if(fin) {
        fin>>n>>f>>g;
        fin.close();
    } else {
        cout<<"Gib n f g ein (z.B. 25 1 0): ";
        cin>>n>>f>>g;
    }

    // 2) System aufbauen
    int dim = (n-1)*(n-1);
    SparseMatrix<double> A = generatePoissonMatrix2D(n);
    VectorXd b              = generateRHS2D(n, f, g);

    // 3) Solver wählen
    cout<<"Direkter LU (d) oder CG (i)? [d/i]: ";
    char c; cin>>c;
    VectorXd x;

    if(c=='d' || c=='D') {
        // 4a) Lösen mit LU
        x = solveSparseLU(A,b);
        // 5a) Speichern als (n-1)×(n-1)
        writeVectorAsMatrix("solution_2_2_1_LU.txt", x, n-1, n-1);
        // 6a) Plot
        //gnuplot_splot("solution_2_2_1_LU.txt");
    }
    else {
        // 4b) Lösen mit CG
        x = solveSparseCG(A,b);
        // 5b) Speichern als (n-1)×(n-1)
        writeVectorAsMatrix("solution_2_2_1_CG.txt", x, n-1, n-1);
        // 6b) Plot
        //gnuplot_splot("solution_2_2_1_CG.txt");
    }

    cout<<"Fertig.\n";
    return 0;
}
