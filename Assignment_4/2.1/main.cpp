#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include "Library.h"
#include "Methods.h"

using namespace std;
using namespace Eigen;

// Function to read matrix data from file
void readMatrix(string filename, MatrixXd &A, VectorXd &b) {
    ifstream file(filename);
    if (!file) {
        cerr << "Error opening file!" << endl;
        exit(1);
    }
    int n;
    file >> n;  // Read matrix size
    A.resize(n, n);
    b.resize(n);
    
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            file >> A(i, j);
        }
    }
    
    for (int i = 0; i < n; i++) {
        file >> b(i);
    }
}

// Function to solve using Dense LU decomposition from the library
void solveDense(MatrixXd &A, VectorXd &b) {
    cout << "Using Dense LU Solver..." << endl;
    VectorXd x = solveDenseLU(A, b);
    writeSingleVector("solution_dense.txt", x);
}

// Function to solve using Sparse LU decomposition from the library
void solveSparse(SparseMatrix<double> &A, VectorXd &b) {
    cout << "Using Sparse LU Solver..." << endl;
    VectorXd x = solveSparseLU(A, b);
    writeSingleVector("solution_sparse.txt", x);
}

int main() {
    string filename = "matrix_data.txt";
    MatrixXd A;
    VectorXd b;
    readMatrix(filename, A, b);
    
    int choice;
    cout << "Select Matrix Type: (1) Dense (2) Sparse: ";
    cin >> choice;
    
    if (choice == 1) {
        solveDense(A, b);
    } else {
        SparseMatrix<double> A_sparse = A.sparseView();
        solveSparse(A_sparse, b);
    }
    return 0;
}
