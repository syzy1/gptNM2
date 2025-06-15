// C++ Standardbibliotheken
#include <string>
#include <iostream>
#include <math.h>
#include <fstream>
#include <cmath>
#include <stdio.h>
#include <iomanip>

// Eigenbibliothek
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/LU>

// Namespaces
using namespace std;
using   Eigen::MatrixXd;
using   Eigen::VectorXd; 

// Bibliothek mit Funktionen
#include "Methods.h"

// LU

#ifndef LIBRARY_H
#define LIBRARY_H

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>
using namespace Eigen;
using namespace std;

// Direct solver for dense
VectorXd solveDenseLU(const MatrixXd &A, const VectorXd &b);

// Direct solver for sparse
VectorXd solveSparseLU(const SparseMatrix<double> &A, const VectorXd &b);

// Iterative solver (e.g. Conjugate Gradient)
VectorXd solveSparseCG(const SparseMatrix<double> &A, const VectorXd &b);

// 2D Poisson matrix (n+1 x n+1 grid => interior is (n-1)*(n-1))
SparseMatrix<double> generatePoissonMatrix2D(int n);

// 2D Poisson right-hand side
VectorXd generateRHS2D(int n, double fVal, double gVal);

#endif

