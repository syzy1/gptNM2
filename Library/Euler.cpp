#include<Library.h>
#include"Euler.h"
#include "Methods.h"

using namespace std;

//Konstruktor
Euler::Euler (double h_euler, double tol_euler, int kmax_euler, void(*function)(VectorXd &, VectorXd &))
{
    // Erstellen der einheitlichen Klasse für alle Euler-Methoden
    //Import der Übergabeparameter; Werden in der Main-Funktion oder einer externen Datei definiert und an den Solver übergeben
    this->h     = h_euler;      //Import der Größe "Schrittweite"
    this->func  = function;     //Import der Funktionen der rechten Seite --> Hier werden Funktionen übergeben, keine konstanten Werte
    this->tol   = tol_euler;    //Import der Toleranzgrenze
    this->kmax  = kmax_euler;   //Import der maximalen Schrittzahl für das numerische Lösungsverfahren
}

//Destructor
Euler::~Euler()
{

}

//Methoden
//Nach Erstellung der Klasse "Euler" mit den oben eingeführten Größen können nun verschiedene Funktionen damit aufgerufen werden.
//Dabei werden weitere Übergabegrößen benötigt. Gemeinsam mit den Standardgrößen sind diese innerhalb der Funktion verwendbar.

//////////////////////////Explizites Eulerverfahren////////////////////////

MatrixXd Euler::explizit(VectorXd y0, VectorXd t_vec)
{
    int n = t_vec.rows();   //Anzahl Zeitschritte aus der Größe des Vektors "t_vec" bestimmen
    int m = y0.rows();      //Anzahl Zustandsgrößen aus der Größe des Vektors "y0" bestimmen
    int i, j;               //Laufvariablen Definieren

    // Speicher Allokierung
    MatrixXd Erg(n, m); 
    VectorXd F(m); 
    VectorXd y_i(m); 

    // Anfangsbedingungen überschreiben
    y_i = y0; 
 
    // Zeitschleife von Startwert bis alle im Zeitvektor "t_vec" defineirten Zeitschritte durchgeführt sind
    for (i=0; i<n; i++)
    {
        // Ergebnisse Speichern 
        for(j=0; j<m; j++)
        {
            Erg(i,j)=y_i(j);    //Schreibt Ergebnisse des aktuellen Zeitschritts in den Gesamt-Lösungsvektor. Dieser wird am Ende der Funktion übergeben. 
        }

        //rechte Seite Lösen
        this->func(y_i, F);     //Lösung der Gleichungen der Rechten Seite mit den aktuellen y-Werten
        //Eulerschritt berechnen
        y_i = y_i + this->h*F;  //Berechnung der Eulerschritts: y(n) = y(n-1) + h*f(y)
      
    }

    //Ergebnismatrix übergeben
    return Erg;  
}

/////////////////////////////////impliziter Euler//////////////////////////////////////
MatrixXd Euler::implizit(VectorXd y0, VectorXd t_vec)
{
    int n   = t_vec.rows();     //Anzahl Zeitschritte aus der Größe des Vektors "t_vec" bestimmen
    int m   = y0.rows();        //Anzahl Zustandsgrößen aus der Größe des Vektors "y0" bestimmen 
    int i,j,k;                  //Laufvariablen definieren
    double h     = this->h;     //Schrittweite von Klassen-Variable auf lokale Variable übertragen
    double h_jac = 1e-8;        //Schrittweite für Jakobimatrix definieren und belegen
    
    //Definition und Allokierung der verwendeten lokalen Matrizen
    MatrixXd Erg((n+1), m); 
    VectorXd F(m); 
    VectorXd y_i(m); 
    VectorXd y_i_plus1(m); 
    VectorXd delta_y(m); 
    MatrixXd Df(y0.size(), y0.size());
    MatrixXd DF(y0.size(), y0.size());
    MatrixXd I(y0.size(), y0.size()); 
    I.setIdentity(); 

    //Anfangsbedingungen überschreiben
    y_i = y0; 
    
    //Zeitschleife
    for(i=0; i<n; i++)
    {
        //Ergebnisse speichern
        for(j=0; j<m; j++)
        {
            Erg(i,j)=y_i(j);    //Schreibt Ergebnisse des aktuellen Zeitschritts in den Gesamt-Lösungsvektor. Dieser wird am Ende der Funktion übergeben. 
        }

        this->func(y_i, F);     //Lösung der Gleichungen der Rechten Seite mit den aktuellen y-Werten

        //expliziter Euler-Schritt
        y_i_plus1 = y_i+h*F; 
        
        k=0;    //Laufvariable für die Newtoniteration innerhalb eines Zeitschritts auf 0 setzten

        //Start der Newton-Iteration im Zeitschritt i
        do{
            this->func(y_i_plus1, F);                   //Auswertung der rechten Seite mit den y-Werten aus dem expliziten Eulerschritt
            Df = jac_cd(this->func, y_i_plus1, h_jac);  //Bildung der Jacobimatrix aus der rechten Seite, den y-Werten aus dem expliziten Eulerschritt und der Schrittweite für die Jakobimatrix (Zentraldifferenz)
            F  = y_i + h*F-y_i_plus1;                   //Berechnung der Matrix F 
            DF = h*Df-I;                                //Berechnung der Matrix DF
            delta_y = DF.fullPivLu().solve(F);          //Lösung des LGS mit LU über Eigenbibliothek 
            y_i_plus1 = y_i_plus1 - delta_y;            //Berechnung der neuen y-Werte durch das erhaltene Delta y

            k++;

        }while (delta_y.norm()>tol*y_i_plus1.norm() && k<kmax);     //Toleranzkriterium prüfen

        //Update neuer Zeitschritt
        y_i=y_i_plus1; 
    }
   
    return Erg; 
}


///////////////////////////SEMI IMPLIZITER EULER /////////////////////////////////////////////

//MatrixXd Euler::semiimplizit(){   ....... hier kommt die Funktion rein ...........}

///////////////////////////SEMI IMPLIZITER EULER mit dense output /////////////////////////////////////////////

///////////////Neville Schema für Richardson Extrapolation /////////////////////////////

//////////////////Steifigkeitsberechnung/////////////////////////////////

// Semi-Implizites Eulerverfahren (Newton-Verfahren in jedem Schritt)
// semiimplizites Euler-Verfahren (implizit + Newton) ohne Extrapolation
MatrixXd Euler::semiimplizit(VectorXd y0, VectorXd t) {
    int n = t.size();
    int m = y0.size();
    MatrixXd Y(n, m);

    // Anfangswert in die erste Zeile schreiben
    Y.row(0) = y0.transpose();

    // Schleife über alle Zeitschritte
    for (int i = 0; i < n - 1; ++i) {
        // Schrittweite aus dem Zeitvektor
        double h_step = t(i+1) - t(i);

        // 1) Expliziter Prädiktor
        VectorXd F(m);
        func(y0, F);
        VectorXd y_new = y0 + h_step * F;

        // 2) Impliziter Korrektur-Schritt via Newton
        for (int iter = 0; iter < kmax; ++iter) {
            VectorXd G(m);
            func(y_new, G);

            // Jacobi-Matrix approximieren
            MatrixXd J = jac_cd(func, y_new, 1e-8);

            // Residuum: y0 + h*G - y_new = 0
            VectorXd res = y0 + h_step * G - y_new;

            // Linearisierte Systemmatrix DF = I - h*J
            MatrixXd DF = MatrixXd::Identity(m, m) - h_step * J;

            // Löse DF * delta = res
            VectorXd delta = DF.lu().solve(res);

            // Aktualisiere y_new
            y_new += delta;

            // Abbruch, wenn Korrektur klein genug
            if (delta.norm() < tol * (1.0 + y_new.norm()))
                break;
        }

        // Für den nächsten Schritt übernehmen
        y0 = y_new;
        Y.row(i + 1) = y0.transpose();
    }

    return Y;
}

// Semi-Implizites Eulerverfahren mit Richardson-Extrapolation
MatrixXd Euler::semiimplizit_do(VectorXd y0, VectorXd t) {
    int n = t.size();
    MatrixXd Y_low = semiimplizit(y0, t);  // Lösung mit Schrittweite h

    double h_backup = h;  // Original h speichern
    h = h / 2.0;
    int n_high = 2 * (n - 1) + 1;
    VectorXd t_high = VectorXd::LinSpaced(n_high, t(0), t(n - 1));
    MatrixXd Y_high = semiimplizit(y0, t_high);  // Lösung mit h/2
    h = h_backup;  // h zurücksetzen

    MatrixXd Y_dense(n, y0.size());
    for (int i = 0; i < n; ++i) {
        // Richardson-Extrapolation (Ordnung erhöhen)
        Y_dense.row(i) = (4.0 * Y_high.row(2 * i) - Y_low.row(i)) / 3.0;
    }

    return Y_dense;
}

// Semi-Implizites Eulerverfahren mit dichter Ausgabe (Interpolation)
MatrixXd Euler::semiimplizit_dense(VectorXd y0,
                                   double tStart,
                                   double tEnd,
                                   const VectorXd &t_dense)
{
    MatrixXd Y_out(t_dense.size(), y0.size());
    int idxDense = 0;

    double t_i = tStart;
    VectorXd y_i = y0;

    // Anfangswerte für t < tStart
    while (idxDense < t_dense.size() && t_dense(idxDense) < tStart) {
        Y_out.row(idxDense) = y0.transpose();
        idxDense++;
    }

    while (t_i < tEnd - 1e-14) {
        double t_ip1 = t_i + h;
        if (t_ip1 > tEnd) t_ip1 = tEnd;
        double h1 = t_ip1 - t_i;

        // Erster Teil-Schritt: y_i -> y_mid
        VectorXd F(y0.size()), G(y0.size());
        func(y_i, F);
        VectorXd y_mid = y_i + 0.5 * h1 * F;

        for (int iter = 0; iter < kmax; iter++) {
            func(y_mid, G);
            VectorXd res = y_i + 0.5 * h1 * G - y_mid;
            MatrixXd J = jac_cd(func, y_mid, 1e-8);
            MatrixXd DF = MatrixXd::Identity(y0.size(), y0.size()) - 0.5 * h1 * J;
            VectorXd delta = DF.lu().solve(res);
            y_mid += delta;
            if (delta.norm() < tol * (1.0 + y_mid.norm())) break;
        }

        // Zweiter Teil-Schritt: y_mid -> y_ip1
        func(y_mid, F);
        VectorXd y_ip1 = y_mid + 0.5 * h1 * F;

        for (int iter = 0; iter < kmax; iter++) {
            func(y_ip1, G);
            VectorXd res = y_mid + 0.5 * h1 * G - y_ip1;
            MatrixXd J = jac_cd(func, y_ip1, 1e-8);
            MatrixXd DF = MatrixXd::Identity(y0.size(), y0.size()) - 0.5 * h1 * J;
            VectorXd delta = DF.lu().solve(res);
            y_ip1 += delta;
            if (delta.norm() < tol * (1.0 + y_ip1.norm())) break;
        }

        // Interpolation mit Parabel durch y_i, y_mid, y_ip1
        while (idxDense < t_dense.size() && t_dense(idxDense) <= t_ip1 + 1e-14) {
            double tau = t_dense(idxDense);
            if (tau < t_i) {
                Y_out.row(idxDense) = y_i;
                idxDense++;
                continue;
            }
            if (tau > t_ip1) break;

            double theta = (tau - t_i) / (t_ip1 - t_i);

            VectorXd y_interp(y0.size());
            for (int comp = 0; comp < y0.size(); comp++) {
                double Yi = y_i(comp);
                double Ymid = y_mid(comp);
                double Yp1 = y_ip1(comp);

                double A = Yp1 - Yi;
                double B = Ymid - Yi;

                double b_ = 4.0 * B - A;
                double a_ = A - b_;
                double c_ = Yi;

                y_interp(comp) = a_ * theta * theta + b_ * theta + c_;
            }
            Y_out.row(idxDense) = y_interp.transpose();
            idxDense++;
        }

        t_i = t_ip1;
        y_i = y_ip1;
    }

    // Für restliche t > tEnd: letzter Wert bleibt konstant
    while (idxDense < t_dense.size()) {
        Y_out.row(idxDense) = y_i.transpose();
        idxDense++;
    }

    return Y_out;
}

// Implementierung von Steifigkeit(...)
double Euler::Steifigkeit(const VectorXd &y0) const {
    // 1) Jacobi-Matrix per zentraler Differenz
    MatrixXd J0 = jac_cd(func, y0, h);

    // 2) alle Eigenwerte (komplex) berechnen
    EigenSolver<MatrixXd> es(J0);
    VectorXcd ev = es.eigenvalues();

    // 3) |Re(λ)| sammeln, Null‐Anteile weglassen
    std::vector<double> absRe;
    for (int i = 0; i < ev.size(); ++i) {
        double rp = ev[i].real();
        if (std::abs(rp) > 1e-12)
            absRe.push_back(std::abs(rp));
    }
    if (absRe.empty()) {
        std::cerr << "[Steifigkeit] alle Eigenwerte ≈ 0!\n";
        return 0.0;
    }

    // 4) SR = max/min
    double vmax = *std::max_element(absRe.begin(), absRe.end());
    double vmin = *std::min_element(absRe.begin(), absRe.end());
    return vmax / vmin;
}