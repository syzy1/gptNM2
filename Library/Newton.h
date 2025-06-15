#include <Library.h>
class Newton
{
public:

    // Konstruktor
    Newton(VectorXd (*)(const VectorXd &), const int = 1000, const double = 1e-10, const double = 1e-5); // <- 1*10^-5 als h, da sonst die Jacobi-Matrix zu ungenau wird (siehe Assignment 0 Aufgabe 5 iii)

    // Methoden
    VectorXd Solve(const VectorXd &);
    MatrixXd JacobiCD1(const VectorXd &, double);

    // Eigenschaften
    int Iterationen;

private:

    // Membervariablen
    VectorXd (*fcn)(const VectorXd &);
    int n_iter_max;
    double tol;
    double h;
};