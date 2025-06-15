#include "Library.h"
#include <iostream>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/Sparse>

using namespace std;
using namespace Eigen;

// ------------------------------------------------
// 1) LU (dense)
VectorXd solveDenseLU(const MatrixXd &A, const VectorXd &b)
{
    if(A.rows()!=A.cols()){
        cerr<<"Matrix not square!\n";
        exit(1);
    }
    return A.fullPivLu().solve(b);
}

// ------------------------------------------------
// 2) Sparse LU
VectorXd solveSparseLU(const SparseMatrix<double> &A, const VectorXd &b)
{
    SparseLU<SparseMatrix<double>> solver;
    solver.analyzePattern(A);
    solver.factorize(A);
    if(solver.info()!=Success){
        cerr<<"Sparse LU factorization failed!\n";
        exit(1);
    }
    return solver.solve(b);
}

// ------------------------------------------------
// 3) Sparse Conjugate Gradient
VectorXd solveSparseCG(const SparseMatrix<double> &A, const VectorXd &b)
{
    ConjugateGradient<SparseMatrix<double>, Lower|Upper> solver;
    solver.compute(A);
    if(solver.info()!=Success){
        cerr<<"CG decomposition failed!\n";
        exit(1);
    }
    return solver.solve(b);
}

// ------------------------------------------------
// 4) Generate 2D Poisson matrix for domain [0,1]x[0,1] with Dirichlet boundary
//    We use n+1 points in each direction => interior has (n-1)*(n-1) unknowns
//    classical 5-point stencil, dimension m=(n-1)*(n-1).
// ------------------------------------------------
SparseMatrix<double> generatePoissonMatrix2D(int n)
{
    // number of interior unknowns
    int m= (n-1)*(n-1);
    vector< Triplet<double> > triplets;
    triplets.reserve(5*m);

    // step size h=1.0/n, factor = 1/h^2
    double h = 1.0/double(n);
    double factor= 1.0/(h*h);

    // We’ll map 2D index (i,j) in [1..n-1]^2 to 1D index: idx= (j-1)*(n-1) + (i-1)
    // where i=1..n-1 is x-direction, j=1..n-1 is y-direction
    // Stencil: 4 on diag, -1 on left/right/up/down
    // Actually that’s: A(k,k)= 4, A(k,k +/-1)= -1, A(k,k +/- (n-1))= -1
    // We multiply that by factor => diag=4*factor, off= -1*factor
    // boundary conditions => no direct addition here, we only do the typical 5-pt interior

    auto idx = [&](int i, int j){
        return (j-1)*(n-1) + (i-1);
    };

    for(int j=1; j<=n-1; j++){
        for(int i=1; i<=n-1; i++){
            int row= idx(i,j);
            // main diag
            triplets.push_back( Triplet<double>(row,row, 4.0*factor) );

            // left neighbor
            if(i>1){
                int col= idx(i-1,j);
                triplets.push_back( Triplet<double>(row,col, -factor) );
            }
            // right neighbor
            if(i< (n-1)){
                int col= idx(i+1,j);
                triplets.push_back( Triplet<double>(row,col, -factor) );
            }
            // down neighbor
            if(j>1){
                int col= idx(i,j-1);
                triplets.push_back( Triplet<double>(row,col, -factor) );
            }
            // up neighbor
            if(j< (n-1)){
                int col= idx(i,j+1);
                triplets.push_back( Triplet<double>(row,col, -factor) );
            }
        }
    }

    SparseMatrix<double> A(m,m);
    A.setFromTriplets(triplets.begin(), triplets.end());
    return A;
}

// ------------------------------------------------
// 5) Generate RHS for 2D PDE: -Laplace u = f, Dirichlet boundary = g
//    We do simplest approach: b= fVal plus terms from boundary (gVal*(some factor)) etc.
//    If fVal=1 => each interior is 1 => multiplied by h^2 => but let's do that explicitly
// ------------------------------------------------
VectorXd generateRHS2D(int n, double fVal, double gVal)
{
    int m= (n-1)*(n-1);
    VectorXd b= VectorXd::Constant(m, fVal);
    // For each interior point (i,j), we might add contributions from boundary if needed:
    // Because boundary is u=gVal => typical approach adds gVal if neighbor is outside domain.
    // We'll do the standard approach: for each (i,j) next to boundary => b(row)+= gVal
    double h= 1.0/double(n);
    double factor= 1.0/(h*h);

    auto idx = [&](int i, int j){
        return (j-1)*(n-1)+(i-1);
    };

    for(int j=1; j<=n-1; j++){
        for(int i=1; i<=n-1; i++){
            int row= idx(i,j);
            // left boundary i=1 => neighbor is i=0 => that’s boundary => + gVal* factor
            if(i==1)   b(row)+= gVal* factor;
            if(i==n-1) b(row)+= gVal* factor;
            if(j==1)   b(row)+= gVal* factor;
            if(j==n-1) b(row)+= gVal* factor;
        }
    }

    return b;
}
