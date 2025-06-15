#include "Library.h"
#include "Methods.h"
#include <fstream>
#include <iostream>

using namespace std;

void writeVectorOneLine(const string &filename, const VectorXd &x)
{
    ofstream out(filename);
    if(!out){
        cerr << "Cannot open output file '"<< filename <<"'\n";
        return;
    }
    // Print all entries on a single line, space-separated
    for(int i=0; i<x.size(); i++){
        out << x(i);
        if(i< x.size()-1) out << " ";
    }
    out << "\n";
    out.close();
    cout << "Wrote solution in one line to " << filename << "\n";
}

// Print x in matrix layout, row by row
void writeVectorAsMatrix(const string &filename,
                         const VectorXd &x,
                         int rows, int cols)
{
    if(rows * cols != x.size()){
        cerr << "[writeVectorAsMatrix] Error: rows*cols != x.size()!\n"
             << "  (rows=" << rows << ", cols=" << cols
             << ", x.size()=" << x.size() << ")\n";
        return;
    }

    ofstream file(filename);
    if(!file.good()){
        cerr << "[writeVectorAsMatrix] Cannot open output file '" << filename << "'!\n";
        return;
    }

    for(int r = 0; r < rows; r++){
        for(int c = 0; c < cols; c++){
            double val = x(r * cols + c);
            file << val;
            // Separate columns by a space (except maybe last column).
            if(c < cols - 1) file << " ";
        }
        file << "\n";  // new line after each row
    }

    file.close();
    cout << "Wrote matrix form to " << filename << "\n";
}

void writeSingleVector(string filename, VectorXd &x)
{
    ofstream file(filename);
    if(!file.good()){
        cerr<<"Cannot open output '"<<filename<<"'\n";
        return;
    }
    for(int i=0; i<x.size(); i++){
        file << x(i) << "\n";
    }
    file.close();
    cout<<"Wrote single-column vector to "<<filename<<"\n";
}

//Tutorial
void   test(int o){


std::cout << "Das Ergebnis ist \t" << o << std::endl;

}

// Zwei Vektoren eine eine Textdatei schreiben fürs Plotten mit gnuplot
void write2Vectors(VectorXd &vec1, VectorXd &vec2, string name)
{
 int i; 
    
    ofstream myfile;
    
    // Datei öffnen
    myfile.open(name);
   
    for (i=0; i<vec1.rows(); i++){
        myfile << vec1(i) << "\t"<< vec2(i)<<endl;
    }

    // Datei schließen
    myfile.close();
}

void write3Vectors(VectorXd &vec1, VectorXd &vec2, VectorXd &vec3, string name)
{
    int i; 
    
    ofstream myfile;
    
    // Datei öffnen
    myfile.open(name);
   
    for (i=0; i<vec1.rows(); i++){
        myfile << vec1(i) << "\t"<< vec2(i)<<"\t"<<vec3(i)<<endl;
    }

    // Datei schließen
    myfile.close();
}

// Vector und Matrix eine eine Textdatei schreiben fürs Plotten mit gnuplot
void writeVectorMatrix(VectorXd &vec, MatrixXd &mat, string name)
{
    //- Datei 
    ofstream myfile;

    myfile.open(name);

    for (int i=0; i<vec.rows(); i++)
    {
        
        myfile << vec(i) << "\t ";

        for (int j = 0; j<mat.cols(); j++)
        {
            myfile << mat(i,j) << "\t ";
        }
        myfile << "\n";
        
    }

    // Datei schließen
    myfile.close(); 
}


// Jacobimatrix Zentraldifferenz
MatrixXd jac_cd(void(*fcn)(VectorXd &, VectorXd &), const VectorXd &x, const double h)
{
    // Deklarationen
    VectorXd xhelp1(x.size());
    VectorXd xhelp2(x.size());
    VectorXd fhelp1(x.size());
    VectorXd fhelp2(x.size());
    MatrixXd jac(x.size(), x.size());

    for (int i=0; i<x.size(); i++) 
    {            
        // Update des x-Vektors
        xhelp1 = x;
        xhelp2 = x;

        xhelp1[i] =  xhelp1[i]-h;
        xhelp2[i] =  xhelp2[i]+h;
                
        //Auswertung der Funktionen an der Stelle x-h
        fcn(xhelp1, fhelp1);
                
        //Auswertung der Funktionen an der Stelle x+h
        fcn(xhelp2, fhelp2);
                
        //Zeile der Jacobi-Matrix df_i/dx_k
        for (int j=0; j<x.size(); j++)
        {     
            jac(j,i) = (fhelp2(j) - fhelp1(j)) / (2.0 * h);     
        }
    }
    return jac;
}


// Frobeniusnorm Matrix
double Frobenius_Mat(const MatrixXd &A)
{
    int i, j;
    double sum = 0;
    double f;

    for (i = 0; i < A.rows(); i++)
    {
        for (j = 0; j < A.cols(); j++)
        {
            sum += A(i, j) * A(i, j); // Summe der Quadrate der einzelnen Einträge
        }
    }

    f = sqrt(sum); // Ergebnis der Frobeniusnorm
    return f;
}

// Frobeniusnorm Vector
double Frobenius_Vec(const VectorXd &b)
{
    int i;
    double sum = 0;
    double f;

    for (i = 0; i < b.size(); i++)
    {
        sum += pow(b(i), 2); // Summe der Quadrate der einzelnen Einträge
    }

    f = sqrt(sum); // Ergebnis der Frobeniusnorm
    return f;
}
