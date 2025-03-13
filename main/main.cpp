#define _USE_MATH_DEFINES
#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <random>
#include <typeinfo>
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
void exp_eigenvals(MatrixXcd eigenvals, double time)
{
  for (int ii=0; ii<eigenvals.rows(); ii++)
  {
    for (int ii=0; ii<eigenvals.rows(); ii++)

  }
}

int main() 
{
  // Definition of parameters 
  vector<int> L = {8, 16, 32, 64}; // # of lattice points
  vector<int> W = {2, 3, 4}; // perturbance amplitude
  const int N = 2; // # of simulations per parameter for the averaging of the simulation values
  const double Tmax = 20; // Total value for the simulation time


  int ii = 0;
  int kk = 0;
  
  // Distribuition function innitialization
  default_random_engine generator;
  uniform_real_distribution<double> distribution(-W[kk], W[kk]);

  // Parameters of the simulation
  vector<int> pars = {L[ii], W[kk], 0};

  // Innitializations of the storage vectors for the hamiltonians and the respective diagonalized matrices
  vector<MatrixXcd> H_vct;
  vector<ComplexEigenSolver<MatrixXcd>> eigen_vctval;
  
  // Innitialization of the initial wave function
  VectorXcd psi_0;
  
  for (int ll=0; ll<N; ll++)
  {
    pars.pop_back();
    pars.push_back(ll);

    // Innitialization of the Matrix object which will be the hamiltonian
    MatrixXcd ham(pars[0], pars[0]);
    
    // Hamiltonian creation with a diagonal potential perturbation
    for (int jj=0; jj<L[ii]-1; jj++)
    {
      ham(jj, jj) = distribution(generator) ;
      ham(jj, jj+1) = 1; 
      ham(jj+1, jj) = 1; 
    }

    H_vct.push_back(ham);

    // Eigen-vector and -values solver
    ComplexEigenSolver<MatrixXcd> ces(ham);

    eigen_vctval.push_back(ces);
    //print_ham_eigen_vectors_values(ham, ces, pars);
    
    // Definition of the initial state of the wave function in the center of the lattice
    psi_0 = VectorXcd::Zero(pars[0])
    psi_0(pars[0]/2) = 1;

    VectorXcd psi(psi_0);
      

    // Initial value of time
    double tt = -5;

    for (double xx=-5; xx<Tmax; xx+=0.5)
    {
      tt += xx*xx;

      psi = ((ces.eigenvector()).conjugate()).transpose() * psi_0;

      
      

      
    }
  }




/*
  Hist/Prob Distr. Test:
  const int nrolls = 100000;  // Number of random samples
  const int nintervals = 1000; // Number of bins in histogram

  std::default_random_engine generator;
  std::gamma_distribution<double> distribution(2.0, 2.0);

  std::vector<double> data; // Store generated random numbers

  for (int i = 0; i < nrolls; ++i) {
      data.push_back(distribution(generator)); // Generate random numbers
  }

  // Create histogram
  plt::hist(data, nintervals, 1);
  plt::title("Uniform Distribution Histogram");
  plt::xlabel("Value");
  plt::ylabel("Frequency");
  plt::show();

  return 0;
*/
}