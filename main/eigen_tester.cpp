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
  
    int n = 4;
    MatrixXcd psi_0(n,n);
    density_creator(psi_0);
    cout << psi_0 << endl << endl;

    default_random_engine generator;
    normal_distribution<double> distribution(2.0,2);
    MatrixXcd ham(n,n);
    for (int i =0; i<1; i++)
    {
        ham.setZero();
        hamiltonian_creator(ham, generator, distribution);
        cout << ham << endl << endl;
        cout << 2*ham/n << endl << endl;
    }
    
}
