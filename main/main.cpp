#define _USE_MATH_DEFINES
#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <random>
#include <complex>
#include <string>
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

void exp_vect_vals(VectorXcd& vect)
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

void mean_vector(vector<vector<double>>& vct, vector<double>& res)
{
  double aux_mean;

  for(int ii=0; ii<vct[0].size(); ii++)
  {
    aux_mean = 0;
    for(int jj=0; jj<vct.size(); jj++)
    {
      aux_mean += vct[jj][ii];
    }

    res.push_back(aux_mean/vct.size());
  }
}

void hamiltonian_creator(MatrixXcd& ham, double dist, bool off_diag=false)
{
  if (off_diag){
    for (int jj=0; jj<ham.rows()-1; jj++)
    {
      ham(jj, jj+1) = dist; 
      ham(jj+1, jj) = dist; 
    }

  } else {
    for (int jj=0; jj<ham.rows()-1; jj++)
    {
      ham(jj, jj) = dist;
      ham(jj, jj+1) = 1; 
      ham(jj+1, jj) = 1; 
    }
    ham(ham.rows()-1,ham.rows()-1) = dist;

  }
}

int main() 
{
  // Definition of parameters 
  vector<int> L = {8, 12, 32, 64}; // # of lattice points
  vector<int> W = {2, 3,4}; // perturbance amplitude
  const int N = 30; // # of simulations per parameter for the averaging of the simulation values
  const double Tmax = 20; // Total value for the simulation time
  bool off_diag = true;

  int kk = 0;
  
  // Distribuition function innitialization
  default_random_engine generator;
  //uniform_real_distribution<double> distribution(-W[kk], W[kk]);
  //gamma_distribution<double> distribution (2.0,2.0);
  normal_distribution<double> distribution (1.0,W[kk]);

  for(int ii = 0; ii<L.size(); ii++)
  {
    cout << "L=" << L[ii] << endl;
    // Parameters of the simulation
    vector<int> pars = {L[ii], W[kk]};

    // Initialization of the initial wave function
    VectorXcd psi_0;

    // Initialization of the vectors where we will store the evolution in time of the value theat characterize every simulated system
    vector<vector<double>> sim_data1;
    vector<vector<double>> sim_data2;

    // Innitialization of the Matrix object which will be the hamiltonian
    MatrixXcd ham;

    for (int ll=0; ll<N; ll++)
    {
      // Hamiltonian creation with a diagonal potential perturbation
      ham.resize(pars[0], pars[0]);
      ham.setZero();
      hamiltonian_creator(ham, distribution(generator));

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

        VectorXcd aux = eigen_vls;
        for(int jjj = 0; jjj < aux.size(); ++jjj)
        {
            aux(jjj) = aux(jjj) * complex<double>(0,-1) * complex<double>(tt,0);
        }
        exp_vect_vals(aux);
        MatrixXcd aux_mat1 = aux.asDiagonal();

        psi = ((ces.eigenvectors()).conjugate()).transpose()*psi_0;
        psi = aux_mat1 * psi;
        psi = ces.eigenvectors() * psi;

        abs2_vect_vals(psi);

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

      sim_data1.push_back(data1);
      sim_data2.push_back(data2);    
    }

    vector<double> graph_data1;
    vector<double> graph_data2;
    vector<double> t_data;

    mean_vector(sim_data1, graph_data1);
    mean_vector(sim_data2, graph_data2);
    for (double xx=-5; xx<Tmax; xx+=0.5)
    {
      t_data.push_back(pow(2, xx));
    }
/*
    for (int oo=0; oo<graph_data1.size(); oo++)
    {
      cout << t_data[oo] << " | " <<graph_data1[oo] << " | " << graph_data2[oo] << endl;
    }
*/
    // Functions plotting
    plt::figure();
    
    plt::suptitle("Normal Distr. Pertur. with L=" + to_string(L[ii]) +" and W=" + to_string(W[kk]));

    // First subplot
    plt::subplot(1, 2, 1); // (rows, columns, index)
    plt::semilogx(t_data, graph_data1, "-x");

    // Second subplot
    plt::subplot(1, 2, 2);
    plt::semilogx(t_data, graph_data2, "-x");

    // Show the figure
    std::string ch1 = to_string(L[ii]);  // Example values
    std::string ch2 = to_string(W[kk]);

    std::string filename = "normal_W=" + ch2 + "_L=" + ch1 + "_off_diag.png";
    plt::save(filename);
    
  }


  return 0;
}