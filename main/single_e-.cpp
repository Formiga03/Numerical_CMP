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
  bool off_diag = false;
  bool period_pot = true;
  
  // Distribution function initialization
  default_random_engine generator;

  for(int kk = 0; kk<W.size(); kk++)
  {
    cout << "W=" << W[kk] << ":" << endl; 

    // Distribution function:
    uniform_real_distribution<double> distribution(0, 2*M_PI);
    //gamma_distribution<double> distribution (2.0,2.0);
    //normal_distribution<double> distribution(2.0,W[kk]);

    // distribuition function and perturbation info
    vector<double> pert = {W[kk]};

    // return probability
    vector<vector<double>> plt_dt1; // data for output in .txt file: {tpars, Lvals, data}
    plt_dt1.push_back(pert);
    plt_dt1.push_back(T_pars);
    plt_dt1.push_back(L);

    // deviation
    vector<vector<double>> plt_dt2; // data for output in .txt file: {tpars, Lvals, data}
    plt_dt2.push_back(pert);
    plt_dt2.push_back(T_pars);
    plt_dt2.push_back(L);

    for(int ii = 0; ii<L.size(); ii++)
    {
      cout << "-L=" << L[ii] << endl;

      // Parameters of the simulation
      vector<int> pars = {int(L[ii]), int(W[kk])};

      // Initialization of the initial wave function
      VectorXcd psi_0;

      // Initialization of the vectors where we will store the evolution in time of the value theat characterize every simulated system
      vector<vector<double>> sim_data1;
      vector<vector<double>> sim_data2;

      // Innitialization of the Matrix object which will be the hamiltonian
      MatrixXcd ham(pars[0], pars[0]);

      for (int ll=0; ll<N; ll++)
      {
        printProgressBar(ll, N);

        hamiltonian_creator(ham, generator, distribution, off_diag);
        
        // Eigen-vector and -values solver
        ComplexEigenSolver<MatrixXcd> ces(ham);
        //print_ham_eigen_vectors_values(ham, ces, pars);
        
        // Definition of the initial state of the wave function in the center of the lattice
        VectorXcd psi_0 = VectorXcd::Zero(pars[0]);
        psi_0(pars[0]/2) = 1;

        VectorXcd psi(psi_0);
        VectorXcd eigen_vls(ces.eigenvalues());
        
        vector<double> data1, data2;

        // Deviation Calc.
        VectorXcd aux_vct1(pars[0]);
        VectorXcd aux_vct2(pars[0]);

        for (int iii = -pars[0]/2; iii < pars[0]/2; ++iii) 
        {
          aux_vct1(iii+pars[0]/2) = complex<double>(iii*iii,0);
          aux_vct2(iii+pars[0]/2) = complex<double>(iii,0);
        }


        // Initial value of time
        double tt;

        for (double xx=T_pars[0]; xx<T_pars[1]; xx+=T_pars[2])
        {
          tt = pow(2, xx);

          VectorXcd aux = eigen_vls;
          for(int jjj = 0; jjj < aux.size(); ++jjj)
          {
              aux(jjj) = aux(jjj) * complex<double>(0,-1) * complex<double>(tt,0);
          }
          exp_vect_vals(aux);
          MatrixXcd aux_mat1 = aux.asDiagonal();

          MatrixXcd U_mat = ces.eigenvectors() * aux_mat1 * ((ces.eigenvectors()).conjugate()).transpose();
          psi =  U_mat * psi_0;

          abs2_vect_vals(psi);

          // Return probability calc.:
          data1.push_back(psi(pars[0]/2).real());
          cout << tt << ":" << endl;
          cout << psi(pars[0]/2).real() << endl;

          // Deviation calc.:
          complex<double> num1 = aux_vct1.dot(psi);
          complex<double> num2 = aux_vct2.dot(psi);

          data2.push_back(num1.real()-num2.real()*num2.real());
          
        }

        sim_data1.push_back(data1);
        sim_data2.push_back(data2);    
      }
      cout<< endl;

      vector<double> graph_data1;
      vector<double> graph_data2;

      mean_vector(sim_data1, graph_data1);
      mean_vector(sim_data2, graph_data2); 
      
      plt_dt1.push_back(graph_data1);
      plt_dt2.push_back(graph_data2);

    }
    cout<< endl;
    iofile writter;
    string filename = "";
    if (not period_pot) {
      filename +=  string(typeid(distribution).name());
    } else {
      filename +=  "QuasiPeriodicPert_withAlpha_";
    }
    if(off_diag){
      filename += "_OffDiag";
    }
    filename += "_W=" + to_string(int(W[kk]));

    writter.output_data("outputs/ReturnProb_"+filename+".txt", plt_dt1);
    writter.output_data("outputs/StandDev_"+filename+".txt", plt_dt2);
  }
  cout<< endl;

  return 0;
}