#define _USE_MATH_DEFINES
#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <random>
#include <complex>
#include <string>
#include <iomanip>
#include <fstream>
#include <chrono>
#include <thread>

#include "../inc/matplotlibcpp.h"
#include "../inc/functions.h"
#include "../inc/iofile.h"

using namespace std;
using namespace Eigen;
namespace plt = matplotlibcpp;

int main() 
{

    Eigen::VectorXd mat(3);
    mat << 1, 2, 3;

    double threshold = 2.0;
    
    // Apply the thresholding operation
    mat = (mat.array() >= pow(10(-20)).select(mat, 0);

    std::cout << "Thresholded Matrix:\n" << mat.array() << std::endl;

  /*
    int n = 4;
    MatrixXcd phi_0(n,n/2);
    VectorXcd aux(n/2);
    for (int jjj=0; jjj<n; jjj++)
    {
        aux.setZero();
        if (jjj%2 == 1){
            aux((jjj-1)/2) = 1;
            cout << jjj << " | " << (jjj+1)/2 << " | " << aux.transpose() << endl ;
            phi_0.row(jjj) = aux;
        }
    }
    cout << phi_0* phi_0.transpose() << endl;

    MatrixXcd psi_0(n,n);
    density_creator(psi_0);

    default_random_engine generator;
    normal_distribution<double> distribution(2.0,2);
    MatrixXcd ham(n,n);
    //ham_PeridPot_creator(ham);
    cout << ham << endl;

    for (int i =0; i<1; i++)
    {
        ham.setZero();
        hamiltonian_creator(ham, generator, distribution);
        cout << ham << endl << endl;
        cout << 2*ham/n << endl << endl;
    }
    */




}
