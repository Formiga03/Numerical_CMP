#define _USE_MATH_DEFINES
#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <random>
#include <complex>
#include "../lib/matplotlibcpp.h"

using namespace std;
using namespace Eigen;
namespace plt = matplotlibcpp;

void eigen_vector(VectorXd vct1, vector<double> vct2)
{
  for (int i = 0; i < vct1.size(); ++i) 
  {
    vct2.push_back(vct1(i));    
  }
}

void print_ham_eigen_vectors_values(MatrixXcd ham, ComplexEigenSolver<MatrixXcd> ces, vector<int> pars)
{
  cout << "Hamiltonian: (L=" << pars[0] << ", W=" << pars[1] << ", Iter=" << pars[2] << ")" << endl;
  cout << ham << endl;
  cout << "_______________________________________________" << endl;
  cout << "Eigenvectors:" << endl;
  cout << ces.eigenvectors() << endl;
  cout << "_______________________________________________" << endl;
  cout << "Eigenvalues:" << endl;
  cout << ces.eigenvalues() << endl;
  cout << "_______________________________________________\n" << endl;

}

void exp_vect_vals(VectorXcd& vect, complex<double> tt)
{
  for(int ii = 0; ii < vect.size(); ++ii)
  {
    vect(ii) = exp(vect(ii));
  }
}

void abs2_vect_vals(VectorXcd& vect)
{
  for(int ii = 0; ii < vect.size(); ++ii)
  {
    vect(ii) = abs(vect(ii))*abs(vect(ii));
  }
}

int main() 
{
const double Tmax = 20; // Total value for the simulation time
vector<int> pars = {4, 3, 0};

MatrixXcd ham(4,4);

ham << std::complex<double>(-0.6093676, 0), std::complex<double>(1.0, 0), std::complex<double>(0.0, 0), std::complex<double>(0.0, 0),
         std::complex<double>(1.0, 0), std::complex<double>(1.58727453, 0), std::complex<double>(1.0, 0), std::complex<double>(0.0, 0),
         std::complex<double>(0.0, 0), std::complex<double>(1.0, 0), std::complex<double>(2.83231254, 0), std::complex<double>(1.0, 0),
         std::complex<double>(0.0, 0), std::complex<double>(0.0, 0), std::complex<double>(1.0, 0), std::complex<double>(-2.77792594, 0);


// Eigen-vector and -values solver
ComplexEigenSolver<MatrixXcd> ces(ham);

//print_ham_eigen_vectors_values(ham, ces, pars);

// Definition of the initial state of the wave function in the center of the lattice
VectorXcd psi_0 = VectorXcd::Zero(pars[0]);
psi_0(pars[0]/2) = 1;

VectorXcd psi(psi_0);
VectorXcd eigen_vls(ces.eigenvalues());

vector<double> data1, data2;

// Initial value of time
double tt;

for (double xx=-5; xx<Tmax; xx+=0.5)
{
    tt = pow(2, xx);
    cout << tt << endl;
    VectorXcd aux = eigen_vls;
    for(int ii = 0; ii < aux.size(); ++ii)
    {
        aux(ii) = aux(ii) * complex<double>(0,-1) * complex<double>(tt,0);
    }
    exp_vect_vals(aux, complex<double>(tt,0));
    MatrixXcd aux_mat1 = aux.asDiagonal();
    //cout << "exp:\n" << aux << endl;
    
    psi = ((ces.eigenvectors()).conjugate()).transpose()*psi_0;
    psi = aux_mat1 * psi;
    psi = ces.eigenvectors() * psi;
    //cout <<"Psi_f:\n" << psi << endl << endl;

    abs2_vect_vals(psi);
    
    //cout << psi << endl;

    data1.push_back(psi(pars[0]/2).real());

    VectorXcd aux_vct1(pars[0]);
    VectorXcd aux_vct2(pars[0]);

    for (int iii = -pars[0]/2; iii < pars[0]/2; ++iii) 
    {
        aux_vct1(iii+pars[0]/2) = complex<double>(iii*iii,0);
        aux_vct2(iii+pars[0]/2) = complex<double>(iii,0);
    }

    complex<double> num1 = aux_vct1.dot(psi);
    complex<double> num2 = aux_vct2.dot(psi);

    data2.push_back(num1.real()-num2.real()*num2.real());
}
for (int oo=0; oo<data1.size(); oo++)
    {
      cout << data1[oo] << " | " << data2[oo] << endl;
    }
}