#include <Library.h>          // enthält Basis-Includes und Eigen
#include "Euler.h"            // semi-implizites Eulerverfahren
#include "Methods.h"          // writeVectorMatrix, usw.
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <string>
#include <cstdlib>
using namespace std;
using namespace Eigen;

// --------------------
// Globale Parameter
// werden direkt in readParameters gesetzt
// --------------------
static double rhoM, rhoP, MW, k;
static double Vplus, Vr0, Vrmax;
static double tStart, tEnd, hStep;
static double Nm0, Np0;

// ---------------------------------------------------------------------
// Liest Parameter aus "rrdat.sec"
// Key:Wert, z.B. "rho_M:800.0"
// ---------------------------------------------------------------------
void readParameters(const string &filename) {
    ifstream fin(filename);
    if (!fin) {
        cerr << "Fehler: Kann '" << filename << "' nicht öffnen\n";
        exit(1);
    }
    string line;
    while (getline(fin, line)) {
        if (line.size()<3) continue;
        string key; double val;
        stringstream ss(line);
        if (getline(ss, key, ':') && (ss >> val)) {
            if      (key=="rho_M")    rhoM   = val;
            else if (key=="rho_P")    rhoP   = val;
            else if (key=="MW")       MW     = val;
            else if (key=="k")        k      = val;
            else if (key=="V_plus")   Vplus  = val;
            else if (key=="V_r_0")    Vr0    = val;
            else if (key=="V_r_max")  Vrmax  = val;
            else if (key=="t_start")  tStart = val;
            else if (key=="t_end")    tEnd   = val;
            else if (key=="h_steps")  hStep  = val;
            else if (key=="NM_0")     Nm0    = val;
            else if (key=="NP_0")     Np0    = val;
        }
    }
    fin.close();
}

// ODE für Modell 1: y=(Nm,Np,Vr)
void polymerOde(VectorXd &y, VectorXd &f) {
    double Nm=y(0), Np=y(1);
    f(0)=Vplus*(rhoM/MW)-k*Nm;
    f(1)=k*Nm;
    f(2)=Vplus + k*MW*Nm*(1.0/rhoP - 1.0/rhoM);
}

// DAE Phase 1: y=(Nm,Np,Vr)
void polymerDae(VectorXd &y, VectorXd &f) {
    double Nm=y(0), Np=y(1), Vr=y(2);
    f(0)=Vplus*(rhoM/MW)-k*Nm;
    f(1)=k*Nm;
    f(2)=-Vr + (MW/rhoM)*Nm + (MW/rhoP)*Np + Vr0;
}

// DAE Phase 2: y=(Nm,Np,Vr), Vplus=0, Vr konstant
void polymerPostSwitch(VectorXd &y, VectorXd &f) {
    double Nm=y(0);
    f(0)=-k*Nm;
    f(1)=k*Nm;
    f(2)=0.0;
}

// Simulation ODE
void simulateOde(const VectorXd &tVec) {
    VectorXd y0(3); y0 << Nm0, Np0, Vr0;
    Euler solver(hStep,1e-6,1000,polymerOde);
    auto Y = solver.semiimplizit_do(y0, tVec);
    writeVectorMatrix(const_cast<VectorXd&>(tVec), Y, "Loesung_A2_3_ODE.txt");
    cout << "ODE in 'Loesung_A2_3_ODE.txt'\n";
}

// Simulation DAE mit Umschaltung
void simulateDaePhaseSwitch(const VectorXd &tVec) {
    int n = tVec.size();
    VectorXd y0(3); y0 << Nm0, Np0, Vr0;
    Euler s1(hStep,1e-6,1000,polymerDae);
    auto Y1 = s1.semiimplizit_do(y0, tVec);
    int idx=-1;
    for(int i=0;i<n;i++) if(Y1(i,2)>=Vrmax){ idx=i; break; }
    if(idx<0) {
        writeVectorMatrix(const_cast<VectorXd&>(tVec), Y1, "Loesung_A2_3_DAE.txt");
        cout<<"Kein Umschalten, gesamte DAE in 'Loesung_A2_3_DAE.txt'\n";
        return;
    }
    VectorXd t1 = tVec.head(idx+1);
    auto Y1s = Y1.topRows(idx+1);
    VectorXd ySw = Y1.row(idx); ySw(2)=Vrmax;
    int n2 = n - idx;
    VectorXd t2(n2);
    for(int i=0;i<n2;i++) t2(i)=tVec(idx+i);
    Euler s2(hStep,1e-6,1000,polymerPostSwitch);
    auto Y2 = s2.semiimplizit_do(ySw, t2);
    for(int i=0;i<Y2.rows();i++) Y2(i,2)=Vrmax;
    int nt = Y1s.rows()+Y2.rows()-1;
    MatrixXd Ytot(nt,3);
    VectorXd tt(nt);
    Ytot.topRows(Y1s.rows()) = Y1s;      tt.head(Y1s.rows()) = t1;
    Ytot.bottomRows(Y2.rows()-1) = Y2.bottomRows(Y2.rows()-1);
    tt.tail(Y2.rows()-1) = t2.tail(Y2.rows()-1);
    writeVectorMatrix(tt, Ytot, "Loesung_A2_3_DAE.txt");
    cout<<"DAE in 'Loesung_A2_3_DAE.txt'\n";
}

int main(){
    readParameters("rrdat.sec");
    int n = int((tEnd-tStart)/hStep)+1;
    VectorXd t = VectorXd::LinSpaced(n, tStart, tEnd);
    simulateOde(t);
    simulateDaePhaseSwitch(t);
    return 0;
}
