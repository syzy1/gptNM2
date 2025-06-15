//Tutorial Funktion
void  test(int); 

// Vektornormen
double Frobenius_Vec(const VectorXd &b);
double Frobenius_Mat(const MatrixXd &A);

// Ausgabefunktionen
void writeVectorMatrix(VectorXd &vec, MatrixXd &mat, string name);
void write2Vectors(VectorXd &vec1, VectorXd &vec2, string name);
void write3Vectors(VectorXd &vec1, VectorXd &vec2, VectorXd &vec3, string name);

// Jakobidifferenz
MatrixXd jac_cd(void(*fcn)(VectorXd &, VectorXd &), const VectorXd &x, const double h);

//LU

#ifndef METHODS_H
#define METHODS_H

#include <Eigen/Dense>
#include <string>
using namespace Eigen;
using namespace std;

// Function to write a single-column vector to a file
void writeSingleVector(string filename, VectorXd &x);

void writeVectorOneLine(const string &filename, const VectorXd &x);
void writeVectorAsMatrix(const std::string &filename,
    const VectorXd &x,
    int rows, int cols);

#endif // METHODS_H


/* evvel beleydi:

#ifndef METHODS_H            // CHANGE: Added include guard
#define METHODS_H            // CHANGE: Added include guard

#include <string>            // CHANGE: Needed for std::string
#include <Eigen/Dense>       // CHANGE: Needed for VectorXd, MatrixXd
using namespace Eigen;       // CHANGE: So we can use VectorXd, MatrixXd directly
using namespace std;         // CHANGE: So we can use 'string' directly

//Tutorial Funktion
void  test(int); 

// Vektornormen
double Frobenius_Vec(const VectorXd &b);
double Frobenius_Mat(const MatrixXd &A);

// Ausgabefunktionen
void writeVectorMatrix(VectorXd &vec, MatrixXd &mat, string name);
void write2Vectors(VectorXd &vec1, VectorXd &vec2, string name);
void write3Vectors(VectorXd &vec1, VectorXd &vec2, VectorXd &vec3, string name);

// Jakobidifferenz
MatrixXd jac_cd(void(*fcn)(VectorXd &, VectorXd &), const VectorXd &x, const double h);

#endif // METHODS_H           // CHANGE: End of include guard

*/