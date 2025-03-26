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
  // Definition of parameters:
  vector<double> L = {8, 12, 32, 64}; // # of lattice points
  vector<double> W = {2, 3, 4}; // perturbance amplitude
  const int N = 30; // # of simulations per parameter for the averaging of the simulation values
  vector<double> T_pars = {-5, 20, 0.5, 1}; // Time specifications of the simulation {tmin, tmax, step, log scale(0:false, 1:trues)}
  bool off_diag = true;
  bool even = true;
  
  // Distribution function initialization
  default_random_engine generator;

  for(int kk = 0; kk<W.size(); kk++)
  {
    cout << "W=" << W[kk] << ":" << endl; 

    // Distribution function:
    //uniform_real_distribution<double> distribution(-W[kk], W[kk]);
    //gamma_distribution<double> distribution (2.0,2.0);
    normal_distribution<double> distribution(2.0,W[kk]);

    // distribuition function and perturbation info
    vector<double> pert = {W[kk]};

    // imbalance
    vector<vector<double>> plt_dt1; // data for output in .txt file: {tpars, Lvals, data}
    plt_dt1.push_back(pert);
    plt_dt1.push_back(T_pars);
    plt_dt1.push_back(L);

    for(int ii = 0; ii<L.size(); ii++)
    {
      cout << "-L=" << L[ii] << endl;

      // Parameters of the simulation
      vector<int> pars = {int(L[ii]), int(W[kk])};

      // Initialization of the initial wave function
      MatrixXcd rho_0(pars[0], pars[0]);

      // Initialization of the density matrux where we will store the evolution in time of the value theat characterize every simulated system
      vector<vector<double>> sim_data1;

      // Innitialization of the Matrix object which will be the hamiltonian
      MatrixXcd ham(pars[0], pars[0]);

      for (int ll=0; ll<N; ll++)
      {
        ham.setZero();
        printProgressBar(ll, N);

        // Hamiltonian creation with a diagonal potential perturbation
        hamiltonian_creator(ham, generator, distribution, off_diag);

        // Eigen-vector and -values solver
        ComplexEigenSolver<MatrixXcd> ces(ham);
        //print_ham_eigen_vectors_values(ham, ces, pars);
        
        // Definition of the initial state of the wave function in the center of the lattice
        density_creator(rho_0, even);

        MatrixXcd rho(rho_0);
        VectorXcd eigen_vls(ces.eigenvalues());
        VectorXcd aux_vct1(pars[0]);
        VectorXcd aux_vct2(pars[0]);
        VectorXcd rho_diag(pars[0]);
        MatrixXcd eigen_vcts(ces.eigenvectors());
        MatrixXcd aux_mat1(pars[0], pars[0]);
        MatrixXcd aux_mat2(pars[0], pars[0]);
        MatrixXcd aux_mat3(pars[0], pars[0]);
        MatrixXcd aux_mat4(pars[0], pars[0]);
        
        vector<double> data1;

        // Initial value of time
        double tt;

        for (double xx=T_pars[0]; xx<T_pars[1]; xx+=T_pars[2])
        {
          tt = pow(2, xx); 
          
          aux_vct1.setZero();
          aux_vct2.setZero();
          rho_diag.setZero();

          aux_mat1.setZero();
          aux_mat2.setZero();
          aux_mat3.setZero();
          aux_mat4.setZero();

          for(int jjj = 0; jjj < pars[0]; ++jjj)
          {
            aux_vct1(jjj) = eigen_vls(jjj) * complex<double>(0,-1) * complex<double>(tt,0);
            aux_vct2(jjj) = eigen_vls(jjj) * complex<double>(0,1) * complex<double>(tt,0);
          }

          exp_vect_vals(aux_vct1);
          exp_vect_vals(aux_vct2);

          aux_mat1 = aux_vct1.asDiagonal(); // diagonalized time evolution exponencial matrix
          aux_mat2 = aux_vct2.asDiagonal(); // hermitian diagonalized time evolution exponencial matrix

          aux_mat3 = eigen_vcts * aux_mat1 * eigen_vcts.conjugate().transpose(); // time evolution mat. 
          aux_mat4 = eigen_vcts * aux_mat2 * eigen_vcts.conjugate().transpose(); // hermitian time evolution mat. 

          rho = aux_mat3 * rho_0 * aux_mat4; 
          rho_diag = rho.diagonal();

          complex<double> even_ind = 0; 
          complex<double> odd_ind = 0; 

          for (int ee = 0; ee < rho_diag.size(); ee += 2) 
          {
            even_ind += rho_diag(ee+1);
            odd_ind += rho_diag(ee);
          }

          data1.push_back(-2*odd_ind.real()/pars[0] + 2*even_ind.real()/pars[0]);

        }

        sim_data1.push_back(data1);
      }
      cout<< endl;

      vector<double> graph_data1;

      mean_vector(sim_data1, graph_data1);
      
      plt_dt1.push_back(graph_data1);

    }
    cout<< endl;
    iofile writter;
    string filename =  string(typeid(distribution).name());
    if(off_diag){
      filename += "_OffDiag";
    }
    filename += "_W=" + to_string(int(W[kk]));

    writter.output_data("outputs/Imbalance_"+filename+".txt", plt_dt1);
  }
  cout<< endl;

  return 0;
}