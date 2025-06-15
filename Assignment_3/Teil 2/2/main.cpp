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
static double Dinf_current;  // hier für C3 (D_{3,∞})

// ODE-Funktion (wie in Teil 2 – Aufgabe 1) für Stoff C3 mit konzentrationsabh. D und Neumann-RB
void ODEConcentrationDiff(VectorXd &y, VectorXd &dydt) {
    int n = y.size();
    dydt = VectorXd::Zero(n);
    for (int i = 0; i < n; ++i) {
        double u_im1 = (i == 0    ? 0.0    : y(i-1));
        double u_i   =              y(i);
        double u_ip1 = (i == n-1 ? y(i) : y(i+1));
        double u_left  = (u_im1 + u_i) / 2.0;
        double u_right = (u_i   + u_ip1) / 2.0;
        double D_left  = Dinf_current + Dinf_current * (1.0 - exp(-0.2 * u_left));
        double D_right = Dinf_current + Dinf_current * (1.0 - exp(-0.2 * u_right));
        double diff_term = (D_right * (u_ip1 - u_i) - D_left * (u_i - u_im1)) / (dx * dx);
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
    // Parameter für Stoff C3
    Dinf_current = 0.0005;
    // Wir simulieren über eine längere Zeit, bis max(u3) < 1
    double t0 = 0.0, tEnd = 20.0;
    double dt = 0.01;
    int nSteps = int((tEnd - t0) / dt) + 1;
    VectorXd tVec(nSteps);
    for (int k = 0; k < nSteps; ++k) tVec(k) = t0 + k * dt;
    int N = int(1.0 / dx) + 1;
    int nInterior = N - 2;

    // Anfangsbedingung für C3: 20 bei x=0.1, sonst 0
    VectorXd y0 = VectorXd::Zero(nInterior);
    int idx_spike = int(0.1 / dx) - 1;
    if (idx_spike >= 0 && idx_spike < nInterior) {
        y0(idx_spike) = 20.0;
    }

    // Lösen des Anfangswertproblems
    Euler solver(dt, 1e-6, 1000, ODEConcentrationDiff);
    MatrixXd Y = solver.semiimplizit(y0, tVec);

    // Maximalwert von u3 zu jedem Zeitschritt ermitteln
    VectorXd max_values(nSteps);
    int drop_index = -1;
    for (int k = 0; k < nSteps; ++k) {
        // Konzentrationsverteilung zu Zeit t_k inkl. Randpunkte rekonstruieren
        double max_u = 0.0;
        // Linker Rand 0, rechter Rand = gleicher Wert wie letzter Innenpunkt (Neumann)
        double u_right = Y(k, nInterior-1);
        if (u_right > max_u) max_u = u_right;
        for (int i = 0; i < nInterior; ++i) {
            double u_val = Y(k, i);
            if (u_val > max_u) max_u = u_val;
        }
        max_values(k) = max_u;
        if (max_u <= 1.0 && drop_index == -1) {
            drop_index = k;
        }
    }

    // Zeitpunkt und Ort, an dem max(u3) erstmals <= 1 fällt:
    if (drop_index != -1) {
        double t_drop = tVec(drop_index);
        // Ort des Maximums zum Zeitpunkt drop_index ermitteln
        // (Suche größtes u3 in der Verteilung)
        double max_val = max_values(drop_index);
        // Alle Positionen der Verteilung zu drop_index betrachten
        int max_pos_index = -1;
        double max_pos_val = -1e9;
        // Randpunkt berücksichtigen
        double u_right = Y(drop_index, nInterior-1);
        if (u_right == max_val) {
            max_pos_index = N-1; // rechter Rand
            max_pos_val = u_right;
        }
        for (int i = 0; i < nInterior; ++i) {
            double u_val = Y(drop_index, i);
            if (u_val == max_val) {
                max_pos_index = i + 1; // +1 wegen Innen->global Index
                max_pos_val = u_val;
                break;
            }
        }
        double x_location = (max_pos_index >= 0 ? max_pos_index * dx : -1.0);
        cout << "Max(u3) unterschreitet 1 bei t = " << t_drop 
             << " (ca.), an Position x = " << x_location << "." << endl;
    } else {
        cout << "Max(u3) ist bis t = " << tEnd << " nicht unter 1 gefallen." << endl;
    }

    // Daten für Plot: max(u3) über der Zeit ausgeben
    write2Vectors(tVec, max_values, string("MaxC3_vs_time.txt"));
}
