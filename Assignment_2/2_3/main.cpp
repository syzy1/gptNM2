#include <Library.h>          // enthält Basis-Includes und Eigen
#include "Euler.h"            // semi-implizites Eulerverfahren
#include "Methods.h"          // writeVectorMatrix, usw.
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <string>
#include <cstdlib>
#include <iomanip>
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
// Liest Parameter aus "rr.dat"
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

// Hilfsfunktion: schreibt Daten mit Kopfzeile in Datei
void writeLabeled(const string &filename, const VectorXd &t, const MatrixXd &Y) {
    ofstream fout(filename);
    fout << "# t\tNm\tNp\tVr\tVdot\n";
    fout << fixed << setprecision(6);
    int n = t.size(), m = Y.cols();
    for (int i = 0; i < n; ++i) {
        fout << t(i);
        for (int j = 0; j < m; ++j) {
            fout << "\t" << Y(i, j);
        }
        fout << "\n";
    }
    fout.close();
}

// ODE für Modell 1 (Befüllung): y=(Nm,Np,Vr)
void polymerOde(VectorXd &y, VectorXd &f) {
    double Nm=y(0), Np=y(1);
    f(0)=Vplus*(rhoM/MW)-k*Nm;
    f(1)=k*Nm;
    f(2)=Vplus + k*MW*Nm*(1.0/rhoP - 1.0/rhoM);
}

// DAE für Modell 2 (Befüllung): y=(Nm,Np,Vr)
void polymerDae(VectorXd &y, VectorXd &f) {
    double Nm=y(0), Np=y(1), Vr=y(2);
    f(0)=Vplus*(rhoM/MW)-k*Nm;
    f(1)=k*Nm;
    f(2)=-Vr + (MW/rhoM)*Nm + (MW/rhoP)*Np + Vr0;
}

// DAE für Modell 1 (Regelung): y=(Nm,Np,Vdot)
void polymerDaeReg1(VectorXd &y, VectorXd &f) {
    double Nm=y(0), Np=y(1), Vd=y(2);
    f(0)=Vd*(rhoM/MW)-k*Nm;                                    // dNm/dt = ...
    f(1)=k*Nm;                                                // dNp/dt = ...
    f(2)=Vd + k*MW*Nm*(1.0/rhoP - 1.0/rhoM);                  // 0 = Vd + ...
}

// ODE für Modell 2 (Regelung, um Index-2 zu umgehen): y=(Nm,Np,Vr)
void polymerOdeReg2(VectorXd &y, VectorXd &f) {
    double Nm=y(0), Np=y(1);
    double Vd = -k*MW*Nm*(1.0/rhoP - 1.0/rhoM);                // Stellgröße
    f(0)=Vd*(rhoM/MW)-k*Nm;                                   // dNm/dt = ...
    f(1)=k*Nm;                                                // dNp/dt = ...
    f(2)=0.0;                                                 // Vr konstant
}

// berechnet durch Parabel-Interpolation den exakten Umschaltzeitpunkt ts,
// sowie Nm(ts), Np(ts)
void computeSwitch(const VectorXd &t, const MatrixXd &Y, int idx,
                   double &ts, double &Nm_s, double &Np_s) {
    int i0 = max(0, idx-2), i1 = max(0, idx-1), i2 = idx;
    double t0 = t(i0), t1 = t(i1), t2 = t(i2);
    double v0 = Y(i0,2), v1 = Y(i1,2), v2 = Y(i2,2);
    Matrix3d M; M << t0*t0, t0, 1,
                   t1*t1, t1, 1,
                   t2*t2, t2, 1;
    Vector3d v(v0, v1, v2);
    Vector3d abc = M.colPivHouseholderQr().solve(v);
    double a=abc(0), b=abc(1), c0=abc(2)-Vrmax;
    double disc = b*b - 4*a*c0;
    double r1 = (-b + sqrt(disc)) / (2*a);
    double r2 = (-b - sqrt(disc)) / (2*a);
    ts = (r1>t0 && r1<t2) ? r1 : r2;
    // Lagrange-Basis für Nm und Np
    auto L = [&](int j, double tx){
        vector<double> tt = {t0,t1,t2};
        double num=1, den=1;
        for(int m=0;m<3;m++) if(m!=j) {
            num *= (tx - tt[m]);
            den *= (tt[j] - tt[m]);
        }
        return num/den;
    };
    Nm_s = Y(i0,0)*L(0,ts) + Y(i1,0)*L(1,ts) + Y(i2,0)*L(2,ts);
    Np_s = Y(i0,1)*L(0,ts) + Y(i1,1)*L(1,ts) + Y(i2,1)*L(2,ts);
}

// Simulation Modell 1: ODE-Befüllung → DAE-Regelung
void simulateModel1(const VectorXd &tVec) {
    int n = tVec.size();
    // 1) Befüllung
    VectorXd y0(3); y0 << Nm0, Np0, Vr0;
    Euler solF(hStep,1e-6,1000,polymerOde);
    MatrixXd Yf = solF.semiimplizit_do(y0, tVec);
    // finde Umschalt-Index
    int idx=-1;
    for(int i=0;i<n;i++) if(Yf(i,2)>=Vrmax){ idx=i; break; }
    // kein Umschalten?
    if(idx<0) {
        MatrixXd Yext(n,4);
        for(int i=0;i<n;i++){
            double Nm_=Yf(i,0), Np_=Yf(i,1), Vr_=Yf(i,2);
            double Vd_ = Vplus;
            Yext(i,0)=Nm_; Yext(i,1)=Np_; Yext(i,2)=Vr_; Yext(i,3)=Vd_;
        }
        writeLabeled("Model1_A2_3.txt", tVec, Yext);
        cout<<"Modell 1 (nur Bef.) in 'Model1_A2_3.txt'\n";
        return;
    }
    // 2) Exakte Umschaltzeit und Werte
    double ts, Nm_s, Np_s;
    computeSwitch(tVec, Yf, idx, ts, Nm_s, Np_s);
    double Vd_s = - k*MW*Nm_s*(1.0/rhoP - 1.0/rhoM);
    // 3) Aufbau Bef.-Ausgabe bis Umschaltpunkt
    MatrixXd Yf_ext(idx+1,4);
    VectorXd t_f(idx+1);
    for(int i=0;i<idx;i++){
        double Nm_=Yf(i,0), Np_=Yf(i,1), Vr_=Yf(i,2);
        double Vd_ = Vplus;
        Yf_ext(i,0)=Nm_; Yf_ext(i,1)=Np_; Yf_ext(i,2)=Vr_; Yf_ext(i,3)=Vd_;
        t_f(i)=tVec(i);
    }
    // Umschalt-Zeile
    Yf_ext(idx,0)=Nm_s; Yf_ext(idx,1)=Np_s;
    Yf_ext(idx,2)=Vrmax; Yf_ext(idx,3)=Vd_s;
    t_f(idx)=ts;
    // 4) Regelung von ts bis Ende
    int m = n - idx;
    VectorXd t_r(m+1);
    t_r(0)=ts;
    for(int j=1;j<=m;j++) t_r(j)=tVec(idx + j -1);
    VectorXd yS(3); yS<<Nm_s,Np_s,Vd_s;
    Euler solR(hStep,1e-6,1000,polymerDaeReg1);
    MatrixXd Yr = solR.semiimplizit_do(yS, t_r);
    // 5) Aufbau Regel-Ausgabe (Vdot im 4.)
    MatrixXd Yr_ext(m+1,4);
    for(int j=0;j<=m;j++){
        double Nm_=Yr(j,0), Np_=Yr(j,1), Vd_=Yr(j,2);
        Yr_ext(j,0)=Nm_; Yr_ext(j,1)=Np_;
        Yr_ext(j,2)=Vrmax; Yr_ext(j,3)=Vd_;
    }
    // 6) Kombinieren und Speichern
    int rows = Yf_ext.rows() + Yr_ext.rows() - 1;
    MatrixXd Ytot(rows,4); VectorXd tt(rows);
    Ytot.topRows(Yf_ext.rows()) = Yf_ext;
    tt.head(Yf_ext.rows()) = t_f;
    Ytot.bottomRows(Yr_ext.rows()-1) = Yr_ext.bottomRows(Yr_ext.rows()-1);
    tt.tail(Yr_ext.rows()-1) = t_r.tail(Yr_ext.rows()-1);
    writeLabeled("Model1_A2_3.txt", tt, Ytot);
    cout<<"Modell 1 (Fill+Reg) in 'Model1_A2_3.txt'\n";
}

// Simulation Modell 2: DAE-Befüllung → ODE-Regelung
void simulateModel2(const VectorXd &tVec) {
    int n = tVec.size();
    // 1) Befüllung
    VectorXd y0(3); y0 << Nm0, Np0, Vr0;
    Euler solF(hStep,1e-6,1000,polymerDae);
    MatrixXd Yf = solF.semiimplizit_do(y0, tVec);
    int idx=-1;
    for(int i=0;i<n;i++) if(Yf(i,2)>=Vrmax){ idx=i; break; }
    if(idx<0) {
        MatrixXd Yext(n,4);
        for(int i=0;i<n;i++){
            double Nm_=Yf(i,0), Np_=Yf(i,1), Vr_=Yf(i,2);
            double Vd_ = Vplus;  // im Befüllbetrieb konstant
            Yext(i,0)=Nm_; Yext(i,1)=Np_; Yext(i,2)=Vr_; Yext(i,3)=Vd_;
        }
        writeLabeled("Model2_A2_3.txt", tVec, Yext);
        cout<<"Modell 2 (nur Bef.) in 'Model2_A2_3.txt'\n";
        return;
    }
    // 2) Exakter Umschaltpunkt
    double ts, Nm_s, Np_s;
    computeSwitch(tVec, Yf, idx, ts, Nm_s, Np_s);
    // 3) Bef.-Ausgabe bis Umschaltpunkt
    MatrixXd Yf_ext(idx+1,4);
    VectorXd t_f(idx+1);
    for(int i=0;i<idx;i++){
        double Nm_=Yf(i,0), Np_=Yf(i,1), Vr_=Yf(i,2);
        Yf_ext(i,0)=Nm_; Yf_ext(i,1)=Np_;
        Yf_ext(i,2)=Vr_; Yf_ext(i,3)=Vplus;
        t_f(i)=tVec(i);
    }
    Yf_ext(idx,0)=Nm_s; Yf_ext(idx,1)=Np_s;
    Yf_ext(idx,2)=Vrmax; Yf_ext(idx,3)=Vplus;
    t_f(idx)=ts;
    // 4) Regelung von ts bis Ende
    int m = n - idx;
    VectorXd t_r(m+1);
    t_r(0)=ts;
    for(int j=1;j<=m;j++) t_r(j)=tVec(idx + j -1);
    double Vd_s = -k*MW*Nm_s*(1.0/rhoP - 1.0/rhoM);
    VectorXd yS(3); yS<<Nm_s,Np_s,Vd_s;
    Euler solR(hStep,1e-6,1000,polymerOdeReg2);
    MatrixXd Yr = solR.semiimplizit_do(yS, t_r);
    // 5) Aufbau Regel-Ausgabe (Vdot aus Zustand berechnen)
    MatrixXd Yr_ext(m+1,4);
    for(int j=0;j<=m;j++){
        double Nm_=Yr(j,0), Np_=Yr(j,1);
        double Vd_ = -k*MW*Nm_*(1.0/rhoP - 1.0/rhoM);
        Yr_ext(j,0)=Nm_; Yr_ext(j,1)=Np_;
        Yr_ext(j,2)=Vrmax; Yr_ext(j,3)=Vd_;
    }
    // 6) Kombinieren und Speichern
    int rows = Yf_ext.rows() + Yr_ext.rows() - 1;
    MatrixXd Ytot(rows,4); VectorXd tt(rows);
    Ytot.topRows(Yf_ext.rows()) = Yf_ext;
    tt.head(Yf_ext.rows()) = t_f;
    Ytot.bottomRows(Yr_ext.rows()-1) = Yr_ext.bottomRows(Yr_ext.rows()-1);
    tt.tail(Yr_ext.rows()-1) = t_r.tail(Yr_ext.rows()-1);
    writeLabeled("Model2_A2_3.txt", tt, Ytot);
    cout<<"Modell 2 (Fill+Reg) in 'Model2_A2_3.txt'\n";
}

int main(){
    readParameters("rr.dat");
    int n = int((tEnd-tStart)/hStep) + 1;
    VectorXd t = VectorXd::LinSpaced(n, tStart, tEnd);

    simulateModel1(t);    // Modell 1: ODE → DAE-Regelung
    simulateModel2(t);    // Modell 2: DAE → ODE-Regelung

    // Erstelle Gnuplot-Skript und rufe Gnuplot auf
    ofstream gp("plot_A2_3.gnuplot");
    gp << R"(# plot_A2_3.gnuplot
# Grafische Auswertung für Aufgabe 2.3
set terminal pngcairo size 1200,900 enhanced font 'Arial,12'
set output 'A2_3_multiplot.png'
set xlabel 'Zeit t [min]'
set ylabel 'Werte'
set grid
set multiplot layout 2,2 title 'Aufgabe 2.3 – Modell 1 vs. Modell 2'

set title 'Monomermenge N_m'
plot \
    'Model1_A2_3.txt' using 1:2 with lines lw 2 title 'Modell 1', \
    'Model2_A2_3.txt' using 1:2 with lines lw 2 title 'Modell 2'

set title 'Polymermenge N_p'
plot \
    'Model1_A2_3.txt' using 1:3 with lines lw 2 title 'Modell 1', \
    'Model2_A2_3.txt' using 1:3 with lines lw 2 title 'Modell 2'

set title 'Reaktorvolumen V_r [m^3]'
plot \
    'Model1_A2_3.txt' using 1:4 with lines lw 2 title 'Modell 1', \
    'Model2_A2_3.txt' using 1:4 with lines lw 2 title 'Modell 2'

set title 'Zulaufstrom \\dot{V}'
plot \
    'Model1_A2_3.txt' using 1:5 with lines lw 2 title 'Modell 1', \
    'Model2_A2_3.txt' using 1:5 with lines lw 2 title 'Modell 2'

unset multiplot
unset output
)";
    gp.close();
    system("gnuplot plot_A2_3.gnuplot");

    return 0;
}
