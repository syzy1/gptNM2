#include<Library.h>
#include<Eigen/Dense>
using namespace Eigen;

class Euler
{
    public:

        // Konstruktor 
        Euler(double h_euler,
              double tol_euler,
              int    kmax_euler,
              void (*function)(VectorXd&, VectorXd&));

        // Destruktor 
        ~Euler();

        // Ihre bestehenden Methoden 
        MatrixXd explizit(VectorXd, VectorXd);
        MatrixXd implizit(VectorXd, VectorXd);
        MatrixXd semiimplizit(VectorXd y0, VectorXd t);
        MatrixXd semiimplizit_do(VectorXd y0, VectorXd t);
        MatrixXd semiimplizit_dense(VectorXd y0,
                                    double   tStart,
                                    double   tEnd,
                                    const VectorXd &t_dense);
        VectorXd Neville();              // nicht implementiert

        // — neu für Aufgabe 2b: Steifigkeitsquotient
        /**
         * Berechnet SR = max|Re(eigen)|/min|Re(eigen)|
         * für J=∂f/∂y bei y0 (Null‐Eigenwerte ignoriert).
         */
        double Steifigkeit(const VectorXd &y0) const;

    private:

        double h;                        // Schrittweite für Jacobi (und evtl. ODE)
        double tol;                      // Toleranz (für Newton etc.)
        int    kmax;                     // max Newton-Iterationen
        VectorXd y0;                     // (wird nicht benutzt von Steifigkeit)
        VectorXd t;                      // (wird nicht benutzt von Steifigkeit)
        void (*func)(VectorXd&,VectorXd&);// RHS-Funktion der DGL
};
