#include <Library.h>

#include "Newton.h"
#include "Methods.h"

Newton::Newton(VectorXd (*Funktion)(const VectorXd &), const int maxIterationsschritte, double Toleranz, double Schrittweite)
{
    this->fcn = Funktion;
    this->n_iter_max = maxIterationsschritte;
    this->tol = Toleranz;
    this->h = Schrittweite;
    this->Iterationen = 0;
};

//Newton-Verfahren mit der L2-Norm des Funktionswertes als Abbruchbedingung

VectorXd Newton::Solve(const VectorXd &x0)
{
    int n = 0; // Iterationsschritte
    VectorXd f_help;
    VectorXd x = x0;
    double res; // Residuum

    f_help = this->fcn(x); // Funtkionsauswertung

    do
    {

        MatrixXd DF = JacobiCD1(x, this->h);
        VectorXd s = DF.colPivHouseholderQr().solve(f_help);
        VectorXd snew = x - s;
        x = snew;
        f_help = this->fcn(x);
        res = f_help.norm();
        n++;
    } while (res > this->tol && n < this->n_iter_max);


        return x;
};

MatrixXd Newton::JacobiCD1(const VectorXd &x, double h)
{

    // Berechnung der Jacobi Matrix
    MatrixXd J(x.size(), x.size());
    VectorXd fx(x.size());
    VectorXd fxph(x.size()); // f(x) plus h
    VectorXd fxmh(x.size()); // f(x) minus h

    // Loop über alle Variablen
    for (int i = 0; i < x.size(); i++)
    {
        // Auswerten der Funktion am aktuellen Punkt
        VectorXd f = this->fcn(x);

        // Erhöhen der i-ten Variable um h
        VectorXd xph = x;
        xph(i) += h;
        VectorXd fxph = this->fcn(xph);

        // Verringern der i-ten Variable um h
        VectorXd xmh = x;
        xmh(i) -= h;
        VectorXd fxmh = this->fcn(xmh);

        // Berechnen der partiellen Ableitung für die i-te Variable
        for (int j = 0; j < x.size(); j++)
        {
            J(j, i) = (fxph(j) - fxmh(j)) / (2 * h);
        }
    }

    return J;
};
