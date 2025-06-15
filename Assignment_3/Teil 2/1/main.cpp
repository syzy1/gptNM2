#include <iostream>
#include <fstream>
#include <cmath>
#include <Eigen/Dense>
#include "Library.h"
#include "Euler.h"
#include "Methods.h"

using namespace std;
using namespace Eigen;

static double v = 0.2;
static double dx = 0.02;
static double Dinf_current;  // aktueller D_{j,∞} für den Stoff

// ODE-Funktion für PDE mit konzentrationsabhängigem D und Neumann-RB am rechten Rand
void ODEConcentrationDiff(VectorXd &y, VectorXd &dydt) {
    int n = y.size();
    dydt = VectorXd::Zero(n);
    for (int i = 0; i < n; ++i) {
        // Nachbarwerte (Neumann-Randbedingung an x=1: u_ip1 = u_i; Dirichlet an x=0: u_im1 = 0)
        double u_im1 = (i == 0    ? 0.0    : y(i-1));
        double u_i   =              y(i);
        double u_ip1 = (i == n-1 ? y(i) : y(i+1));
        // Lokale Diffusionskoeffizienten an Zellgrenzen (abhängig vom Konzentrationsmittel)
        double u_left  = (u_im1 + u_i) / 2.0;
        double u_right = (u_i   + u_ip1) / 2.0;
        double D_left  = Dinf_current + Dinf_current * (1.0 - exp(-0.2 * u_left));
        double D_right = Dinf_current + Dinf_current * (1.0 - exp(-0.2 * u_right));
        // Diffusionsfluss-Differenz: (D_right*(u_ip1 - u_i) - D_left*(u_i - u_im1)) / (Δx^2)
        double diff_term = (D_right * (u_ip1 - u_i) - D_left * (u_i - u_im1)) / (dx * dx);
        // Konvektionsterm (Upwind-Schema, v > 0 -> Rückwärtsdifferenz)
        double conv_term;
        if (v > 0) {
            conv_term = v * (u_i - u_im1) / dx;
        } else {
            conv_term = v * (u_ip1 - u_i) / dx;
        }
        dydt(i) = diff_term - conv_term;
    }
}

int main() {
    // Zeit- und Raumdiskretisierung
    double t0 = 0.0, tEnd = 1.0;
    double dt = 0.01;
    int nSteps = int((tEnd - t0) / dt) + 1;
    VectorXd tVec(nSteps);
    for (int k = 0; k < nSteps; ++k) tVec(k) = t0 + k * dt;
    int N = int(1.0 / dx) + 1;
    int nInterior = N - 2;
    VectorXd x_all = VectorXd::LinSpaced(N, 0.0, 1.0);

    // D_{j,∞}-Werte für die drei Stoffe (gegeben)
    double Dinf[4] = {0, 0.01, 0.005, 0.0005};

    // Simulation für Stoff 1, 2, 3
    for (int j = 1; j <= 3; ++j) {
        Dinf_current = Dinf[j];
        // Anfangsbedingung: 20 an x=0.1, sonst 0
        VectorXd y0 = VectorXd::Zero(nInterior);
        int idx_spike = int(0.1 / dx) - 1;
        if (idx_spike >= 0 && idx_spike < nInterior) {
            y0(idx_spike) = 20.0;
        }

        Euler solver(dt, 1e-6, 1000, ODEConcentrationDiff);
        MatrixXd Y = solver.semiimplizit(y0, tVec);
        // Konzentrationsprofil am Endzeitpunkt t = 1.0 (mit Randpunkten)
        VectorXd u_end = VectorXd::Zero(N);
        u_end(0) = 0.0;                        // Dirichlet-Bedingung
        u_end(N-1) = Y(nSteps-1, nInterior-1); // Neumann-Bedingung: u(N-1) = u(N-2)
        for (int i = 1; i < N-1; ++i) {
            u_end(i) = Y(nSteps-1, i-1);
        }

        // Ergebnis in Datei schreiben (x und Konzentration bei t=1)
        string filename = "Stoff" + to_string(j) + "_t1.txt";
        MatrixXd out(N, 1);
        out.col(0) = u_end;
        writeVectorMatrix(x_all, out, filename);
    }
}
