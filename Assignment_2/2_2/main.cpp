#include <Library.h>    // Basis-Includes und Eigen
#include "Euler.h"      // semi-implizites Eulerverfahren
#include "Methods.h"    // writeVectorMatrix
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cstdlib>

using namespace std;
using namespace Eigen;

// --------------------
// Globale Parameter
// werden in readParameters gesetzt
// --------------------
static int    modelChoice;
static double rhoM, rhoP, MW, k;
static double Vr0, Nm0, Np0;
static double tStart, tEnd, hStep;
static double Vplus;

// ---------------------------------------------------------------------
// Liest Parameter aus "rr.dat"
// Format: Key:Wert, z.B. "rho_M:800.0"
// ---------------------------------------------------------------------
void readParameters(const string &fname) {
    ifstream fin(fname);
    if (!fin) {
        cerr << "Fehler: Kann '" << fname << "' nicht öffnen\n";
        exit(1);
    }
    string line;
    while (getline(fin, line)) {
        if (line.size()<3) continue;
        string key; double val;
        stringstream ss(line);
        if (getline(ss, key, ':') && (ss >> val)) {
            if      (key=="rho_M")   rhoM     = val;
            else if (key=="rho_P")   rhoP     = val;
            else if (key=="MW")      MW       = val;
            else if (key=="k")       k        = val;
            else if (key=="V_r_0")   Vr0      = val;
            else if (key=="NM_0")    Nm0      = val;
            else if (key=="NP_0")    Np0      = val;
            else if (key=="t_start") tStart   = val;
            else if (key=="t_end")   tEnd     = val;
            else if (key=="h_steps") hStep    = val;
            else if (key=="V_plus")  Vplus    = val;
        }
    }
    fin.close();
}

// ---------------------------------------------------------------------
// Liest Messzeiten aus "exp.txt" und filtert 0 < t ≤ tEnd
// ---------------------------------------------------------------------
VectorXd readExpTimes(const string &filename) {
    ifstream fin(filename);
    if (!fin) {
        cerr << "Fehler: Kann '" << filename << "' nicht öffnen\n";
        exit(1);
    }
    vector<double> T;
    string line;
    while (getline(fin, line)) {
        if (line.empty()) continue;
        double t = stod(line);
        if (t > tStart && t <= tEnd) T.push_back(t);
    }
    fin.close();
    sort(T.begin(), T.end());
    VectorXd out(T.size());
    for (size_t i=0; i<T.size(); ++i) out(i) = T[i];
    return out;
}

// ---------------------------------------------------------------------
// ODE-System (Modell 1): y = (Nm, Np, Vr)
// ---------------------------------------------------------------------
void polymerOde(VectorXd &y, VectorXd &f) {
    double Nm=y(0), Np=y(1);
    f = VectorXd::Zero(3);
    f(0) = Vplus*(rhoM/MW) - k*Nm;
    f(1) = k*Nm;
    f(2) = Vplus + k*MW*Nm*(1.0/rhoP - 1.0/rhoM);
}

// ---------------------------------------------------------------------
// DAE-System (Modell 2): y = (Nm, Np, Vr)
// ---------------------------------------------------------------------
void polymerDae(VectorXd &y, VectorXd &f) {
    double Nm=y(0), Np=y(1), Vr=y(2);
    f = VectorXd::Zero(3);
    f(0) = Vplus*(rhoM/MW) - k*Nm;
    f(1) = k*Nm;
    f(2) = -Vr + (MW/rhoM)*Nm + (MW/rhoP)*Np + Vr0;
}

int main() {
    // 1) Parameter aus Datei einlesen
    readParameters("rr.dat");

    // 2) Messzeiten einlesen & prüfen
    VectorXd tDense = readExpTimes("exp.txt");
    if (tDense.size() == 0) {
        cerr << "Keine Messzeiten gefunden. Abbruch.\n";
        return 1;
    }

    // 3) Modell wählen
    cout << "Modell auswaehlen (1=ODE, 2=DAE): ";
    cin >> modelChoice;
    if (!cin || (modelChoice!=1 && modelChoice!=2)) {
        cout << "Ungueltig, verwende ODE (1).\n";
        modelChoice = 1;
    }

    // 4) Startwerte setzen
    VectorXd y0(3);
    y0 << Nm0, Np0, Vr0;

    // 5) Solver konfigurieren
    double tol = 1e-6;
    int    kmax = 1000;
    auto rhs = (modelChoice==1 ? polymerOde : polymerDae);
    Euler solver(hStep, tol, kmax, rhs);

    // 6) Dichte Lösung berechnen
    MatrixXd Y = solver.semiimplizit_dense(y0, tStart, tEnd, tDense);

    // 7) Ergebnis speichern mit Header
    string outFile = (modelChoice==1
                      ? "Loesung_A2_2_ODE.txt"
                      : "Loesung_A2_2_DAE.txt");
    ofstream out(outFile);
    // Header-Zeile
    out << "t\tNm\tNp\tVr\n";
    // Daten schreiben
    for (int i = 0; i < Y.rows(); ++i) {
        out << tDense(i) << "\t"
            << Y(i,0)    << "\t"
            << Y(i,1)    << "\t"
            << Y(i,2)    << "\n";
    }
    out.close();

    cout << "Ergebnis in '" << outFile << "' gespeichert ("
         << tDense.size() << " Punkte).\n";
    return 0;
}
