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
#include <omp.h> 

#include "../inc/functions.h"
#include "../inc/iofile.h"

using namespace std;
using namespace Eigen;

int main() 
{
  
  // Definition of parameters:
  cout << "Insert the system's lattice sizes you wish to analyse:\n";
  vector<double> L = readVectorFromInput(); // # of lattice points
  cout << "Insert the system's perturbation amplitude you wish to analyse:\n";
  vector<double> W = readVectorFromInput(); // perturbation amplitude
  cout << "Random potencial field in each lattice (false) or random hopping potential (true):\n";
  bool off_diag;
  cin >> off_diag; // Parameter which controls if the hamiltonian has an random potential
                   // field at each lattice point (false) or if the hopping potential is
                   // randomized (true)
  bool PBC=false;
  const int N = 30; // # of simulations per parameter for the averaging of the simulation 
                    // values
  vector<double> T_pars = {-5, 20, 0.5, 1}; // Time specifications of the simulation
                                            // {tmin, tmax, step, log scale(0:false, 1:trues)}
    
  // Distribution function initialization
  default_random_engine generator;

  // data for output in .txt file: {Wval, tpars, Lvals, data}
  vector<vector<double>> plt_dt1, plt_dt2, plt_dt3;
  
  for(int kk = 0; kk<W.size(); kk++)
  {
    cout << "W=" << W[kk] << ":" << endl; 

    // Distribution function:
    uniform_real_distribution<double> distribution(0, 2*M_PI);

    // distribuition function and perturbation info
    vector<double> pert = {W[kk]};

    // return probability
    plt_dt1.push_back(pert);
    plt_dt1.push_back(T_pars);
    plt_dt1.push_back(L);

    // standard deviation
    plt_dt2.push_back(pert);
    plt_dt2.push_back(T_pars);
    plt_dt2.push_back(L);

    // IPR:
    plt_dt3.push_back(pert);
    plt_dt3.push_back(T_pars);
    plt_dt3.push_back(L);
    

    for(int ii = 0; ii<L.size(); ii++)
    {
      cout << "-L=" << L[ii] << endl;
      printProgressBar(ii, L.size());

      // Parameters of the simulation
      vector<int> pars = {int(L[ii]), int(W[kk])};

      // Initialization of the initial wave function
      VectorXcd psi_0;

      // Initialization of the vectors where we will store the evolution in time of the value theat characterize every simulated system
      vector<vector<double>> sim_data1, sim_data2, sim_data3;

      // Innitialization of the Matrix object which will be the hamiltonian
      MatrixXcd ham(pars[0], pars[0]);

      //#pragma omp parallel for
      for (int ll=0; ll<N; ll++)
      {
        // Hamiltonian creation with a diagonal potential perturbation
        ham_PeridPot_creator(ham, W[kk], distribution(generator), off_diag, PBC);

        // Eigen-vector and -values solver
        ComplexEigenSolver<MatrixXcd> ces(ham);
        
        // Definition of the initial state of the wave function in the center of the lattice
        psi_0 = VectorXcd::Zero(pars[0]);
        psi_0(pars[0]/2) = 1;

        // Initialization of auxiliary vectors
        VectorXcd psi(pars[0]);
        VectorXcd eigen_vls(ces.eigenvalues());
        VectorXcd aux_vct1(pars[0]), aux_vct2(pars[0]);

        for (int iii = -pars[0]/2; iii < pars[0]/2; ++iii) 
        {
        aux_vct1(iii+pars[0]/2) = complex<double>(iii*iii,0);
        aux_vct2(iii+pars[0]/2) = complex<double>(iii,0);
        }
        
        vector<double> data1, data2, data3;

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

          // Return probability calc.:
          abs2_vect_vals(psi);
          data1.push_back(psi(pars[0]/2).real());

          //Standad Deviation calc.:
          complex<double> num1 = aux_vct1.dot(psi);
          complex<double> num2 = aux_vct2.dot(psi);
          data2.push_back(num1.real()-num2.real()*num2.real()); 
          
          //IPR calc.:
          data3.push_back((psi.asDiagonal()*psi).sum().real());
        }
        // Adding calculated data from the time progression into vectors to be averaged 
        // over the number of different hamiltonians created
        sim_data1.push_back(data1);
        sim_data2.push_back(data2); 
        sim_data3.push_back(data3); 
      }
      cout<< endl;
      // Vectors which will store the averaged values for different lattice size values for a
      // certain perturbation value
      vector<double> graph_data1, graph_data2, graph_data3;

      // Averaging for the number of hamiltonians with different random phase
      mean_vector(sim_data1, graph_data1);
      mean_vector(sim_data2, graph_data2); 
      mean_vector(sim_data3, graph_data3); 
      
      plt_dt1.push_back(graph_data1);
      plt_dt2.push_back(graph_data2);
      plt_dt3.push_back(graph_data3);

    }
    cout<< endl;

    // Creating .txt file where the simulation's information will be stored to be plot later
    // with the main.py. The text file name has the following structure:
    // "QuantatyName_QuasiPeriodicPert_withAlpha_NatureOfThePotencial(diagonal or off diagonal)_
    //  IfPBCorNot_W=PerturbationValue.txt"
    iofile writter;
    string filename = "";
    filename +=  "QuasiPeriodicPert_withAlpha_";
    
    if(off_diag){
      filename += "_OffDiag";
    }

    if(off_diag){
      filename += "_PBC";
    }

    writter.output_data("outputs/ReturnProb_"+filename+"_W=" + to_string(int(W[kk]))+".txt", plt_dt1);
    writter.output_data("outputs/StandDev_"+filename+"_W=" + to_string(int(W[kk]))+".txt", plt_dt2);
    writter.output_data("outputs/IPR_"+filename+"_W=" + to_string(int(W[kk]))+".txt", plt_dt3);
    
    plt_dt1.clear();
    plt_dt2.clear();
    plt_dt3.clear();
  }

  cout<< endl;

  return 0;
}