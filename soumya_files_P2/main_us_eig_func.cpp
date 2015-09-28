#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;
#include <iomanip>
using std::setw;
#include "armadillo"
using namespace arma;
#include <time.h>

ofstream ofile;

int main(int argc, char *argv[])
{
    double rho_min = 0;
    double rho_max = 5;
    double epsilon;
    epsilon = 1.0e-8;
    int n;
    cout << "Please enter value of n:\n>";
    cin >> n;
    cout << "n = " << n << endl;

    double h = (rho_max - rho_min)/(n+1); //step length
    double w = 0.01;

    vec V(n);
    for (int i = 0; i < n; i++)
    {
        V(i) = pow(rho_min + (i+1)*h,2);
    }

    mat A(n,n);
    A.zeros();

    for (int i=0; i<n; i++)
    {
        A(i,i) = 2/pow(h,2) + V(i);
    }

    double off_diagonal;
    off_diagonal = -1/pow(h,2);

    for (int i=1; i<n; i++)
    {
        A(i,i-1) = off_diagonal;
    }
    for (int i=0; i<n-1; i++)
    {
        A(i,i+1) = off_diagonal;
    }

     //cout << "A: "<< A <<endl;



    vec eigval;
    mat eigvec;
    clock_t start, finish;
    start = clock();
    eig_sym(eigval, eigvec, A);
    finish = clock();
    //double clocks_per_sec
    double t = ((finish-start)/CLOCKS_PER_SEC);

    cout << eigval(0)<<endl;
    cout << eigval(1)<<endl;
    cout << eigval(2)<<endl;
    cout << CLOCKS_PER_SEC;
    cout<< "start: " << start << endl;
    cout << "elapsed time: "<< t <<"s"<< endl;
    cout << "finish: "<< finish <<endl;
    //cout << eigvec << endl;
    //Printing results to file
        ofstream myfile ("t_for_eigfunc.txt"); //Creates output file
          if (myfile.is_open())           //checkes whether the output file is open.
                                          //if open, the following things are put in the output file
          {
              myfile << "n = " << n << endl;
              myfile << "rho_max = "<< rho_max << endl;
              myfile << "min1 = "<< eigval(0) << endl;
              myfile << "min2 = "<< eigval(1) << endl;
              myfile << "min3 = "<< eigval(2) << endl;
              myfile << "elapsed time: = "<< t <<"s"<< endl;
          myfile.close();
          }
          else cout << "Unable to open file";

    return 0;
}
