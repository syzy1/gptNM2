#include <Library.h>       // Basis-Includes und Eigen
#include "Euler.h"         // semi-implizites Eulerverfahren
#include "Methods.h"       // writeVectorMatrix
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>           // für system()
using namespace std;
using namespace Eigen;

// --------------------
// Globale Parameter
// werden alle in readParameters gesetzt
// --------------------
static int    modelChoice;
static double rhoM, rhoP, MW, k;
static double Vr0, Vrmax, Nm0, Np0;
static double tStart, tEnd, hStep;
static double Vplus;

// Liest alle Parameter aus "rrdat.sec"
// Format jeder Zeile: Schlüssel:Wert
void readParameters(const string &file) {
    ifstream fin(file);
    if (!fin) {
        cerr << "Fehler: Kann '" << file << "' nicht öffnen\n";
        exit(1);
    }
    string line;
    while (getline(fin, line)) {
        if (line.size() < 3) continue;
        string key; double val;
        stringstream ss(line);
        if (getline(ss, key, ':') && (ss >> val)) {
            if      (key=="rho_M")    rhoM    = val;
            else if (key=="rho_P")    rhoP    = val;
            else if (key=="MW")       MW      = val;
            else if (key=="k")        k       = val;
            else if (key=="V_plus")   Vplus   = val;
            else if (key=="V_r_0")    Vr0     = val;
            else if (key=="V_r_max")  Vrmax   = val;
            else if (key=="NM_0")     Nm0     = val;
            else if (key=="NP_0")     Np0     = val;
            else if (key=="t_start")  tStart  = val;
            else if (key=="t_end")    tEnd    = val;
            else if (key=="h_steps")  hStep   = val;
        }
    }
    fin.close();
}

// ODE-System (Modell 1): y = (Nm, Np, Vr)
void polymerOde(VectorXd &y, VectorXd &f) {
    double Nm = y(0), Np = y(1);
    f(0) = Vplus*(rhoM/MW) - k*Nm;
    f(1) = k*Nm;
    f(2) = Vplus + k*MW*Nm*(1.0/rhoP - 1.0/rhoM);
}

// DAE-System (Modell 2): y = (Nm, Np, Vr)
void polymerDae(VectorXd &y, VectorXd &f) {
    double Nm = y(0), Np = y(1), Vr = y(2);
    f(0) = Vplus*(rhoM/MW) - k*Nm;
    f(1) = k*Nm;
    f(2) = -Vr + (MW/rhoM)*Nm + (MW/rhoP)*Np + Vr0;
}

int main() {
    // 1) Parameter laden
    readParameters("rrdat.sec");

    // 2) Modell wählen
    cout << "Modell waehlen: 1=ODE, 2=DAE: ";
    cin >> modelChoice;
    if (modelChoice!=1 && modelChoice!=2) {
        cout << "Ungueltig, nehme ODE (1)\n";
        modelChoice = 1;
    }

    // 3) Zeitgitter anlegen
    int nSteps = int((tEnd - tStart)/hStep) + 1;
    VectorXd tVec = VectorXd::LinSpaced(nSteps, tStart, tEnd);

    // 4) Anfangs­werte
    VectorXd y0(3);
    y0 << Nm0, Np0, Vr0;
    cout << "Start: Nm="<<Nm0<<", Np="<<Np0<<", Vr="<<Vr0<<"\n";

    // 5) Steifigkeit prüfen
    auto rhs = (modelChoice==1 ? polymerOde : polymerDae);
    Euler tester(hStep,1e-6,1000,rhs);
    double SR = tester.Steifigkeit(y0);
    cout << "Steifigkeitsquotient SR="<<SR<<"\n";
    if (SR>1e3) cout<<"Hinweis: System ist steif\n";

    // 6) Rechnung mit semi-implizitem Euler
    Euler solver(hStep,1e-6,1000,rhs);
    MatrixXd Y = solver.semiimplizit_do(y0, tVec);

    // 7) Ergebnis speichern
    string out = (modelChoice==1
                  ? "Loesung_A2_1_ODE.txt"
                  : "Loesung_A2_1_DAE.txt");
    writeVectorMatrix(tVec, Y, out);
    cout<<"Fertig. Datei: "<<out<<"\n";
    cout<<"Endwerte: Nm="<<Y(nSteps-1,0)
        <<", Np="<<Y(nSteps-1,1)
        <<", Vr="<<Y(nSteps-1,2)<<"\n";

    // 8) Gnuplot-Skript generieren und ausführen
    string gpFile = (modelChoice==1
                     ? "plot_A2_1_ODE.gp"
                     : "plot_A2_1_DAE.gp");
    ofstream gp(gpFile);
    gp << "reset\n"
       << "set title 'Loesung " << (modelChoice==1 ? "ODE" : "DAE")
          << ": Semi-implizites Euler + Extrapolation'\n"
       << "set xlabel 't [min]'\n"
       << "set ylabel 'Nm, Np, Vr'\n"
       << "plot '" << out << "' using 1:2 with lines title 'Nm', \\" << "\n"
       << "     '" << out << "' using 1:3 with lines title 'Np', \\" << "\n"
       << "     '" << out << "' using 1:4 with lines title 'Vr'\n";
    gp.close();
    system((string("gnuplot ")+gpFile+" -").c_str());

    return 0;
}
