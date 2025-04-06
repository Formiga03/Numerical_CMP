#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <sstream>
#include <Eigen/Dense>
#include <vector>
#include <complex>
#include <string>
#include <iomanip>
#include <random>
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
  cout << "System: Unidimensional lattice with L size and W perturbation amplitude with a single eletron in the middle:\n";
  cout << "Insert the system's lattice sizes you wish to analyse:\n";
  vector<double> L = readVectorFromInput(); // # of lattice points
  cout << "Insert the system's perturbation amplitude you wish to analyse:\n";
  vector<double> W = readVectorFromInput(); // perturbation amplitude
  
  const int N = 30; // # of simulations per parameter for the averaging of the simulation values
  vector<double> T_pars = {-5, 20, 0.5, 1}; // Time specifications of the simulation {tmin, tmax, step, log scale(0:false, 1:trues)}
  cout << "Random potencial field in each lattice (false) or random hopping potential (true):\n";
  bool off_diag;
  cin >> off_diag;
  
  // Distribution function initialization
  default_random_engine generator;

  #pragma omp parallel for
  for(int kk = 0; kk<W.size(); kk++)
  {
    // Distribution function:
    //uniform_real_distribution<double> distribution(0, 2*M_PI);
    normal_distribution<double> distribution(2.0,W[kk]);

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
      // Parameters of the simulation
      vector<int> pars = {int(L[ii]), int(W[kk])};

      // Initialization of the vectors where we will store the evolution in time of the value
      // that characterizes every simulated system.
      vector<vector<double>> sim_data1, sim_data2;

      // Initialization of the Matrix object which will be the hamiltonian
      MatrixXcd ham(pars[0], pars[0]);
      ham.setZero();

      //Initialization of auxiliary vectors
      VectorXcd psi_0(pars[0]);
      VectorXcd aux_vct1(pars[0]), aux_vct2(pars[0]);
      VectorXcd psi(pars[0]), eigen_vls(pars[0]);

      // Auxliary vector for the Standard Deviation Calc.:
      for (int iii = -pars[0]/2; iii < pars[0]/2; ++iii) 
      {
        aux_vct1(iii+pars[0]/2) = complex<double>(iii*iii,0);
        aux_vct2(iii+pars[0]/2) = complex<double>(iii,0);
      }
      
      for (int ll=0; ll<30; ll++)
      {
        printProgressBar(ll, N);
        //multiLevelProgressBar(kk, W[kk], ii, L[ii], ll, N);

        hamiltonian_creator(ham, generator, distribution, off_diag);
        
        // Eigen-vector and -values solver
        ComplexEigenSolver<MatrixXcd> ces(ham);
        
        // Definition of the initial state of the wave function in the center of the lattice
        psi_0.setZero();
        psi.setZero();
        psi_0(pars[0] / 2) = 1; 
        eigen_vls = ces.eigenvalues();
        
        vector<double> data1, data2;
        
        for (double xx=T_pars[0]; xx<T_pars[1]; xx+=T_pars[2])
        {
          double tt = exp2(xx);

          VectorXcd aux = eigen_vls;
          for(int jjj = 0; jjj < aux.size(); ++jjj)
          {
              aux(jjj) = aux(jjj) * complex<double>(0,-1) * complex<double>(tt,0);
          }
          aux = aux.array().exp();
          MatrixXcd aux_mat1 = aux.asDiagonal();

          MatrixXcd U_mat = ces.eigenvectors() * aux_mat1 * ((ces.eigenvectors()).conjugate()).transpose();
          psi =  U_mat * psi_0;

          // Return probability calc.:
          psi = psi.array().cwiseAbs().pow(2);
          data1.push_back(psi(pars[0]/2).real());

          //Standad Deviation calc.:
          complex<double> num1 = aux_vct1.dot(psi);
          complex<double> num2 = aux_vct2.dot(psi);

          data2.push_back(num1.real()-num2.real()*num2.real()); 
          
        }
        // Adding calculated data from the time progression into vectors to be averaged 
        // over the number of different hamiltonians created
        sim_data1.push_back(data1);
        sim_data2.push_back(data2);    
      }
      // Vectors which will store the averaged values for different lattice size values for a
      // certain perturbation value
      vector<double> graph_data1, graph_data2;

      // Averaging for the number of hamiltonians with different random potencial fields
      mean_vector(sim_data1, graph_data1);
      mean_vector(sim_data2, graph_data2); 
      
      plt_dt1.push_back(graph_data1);
      plt_dt2.push_back(graph_data2);
    }

    // Creating .txt file where the simulation's information will be stored to be plot later
    // with the main.py. The text file name has the following structure:
    // "QuantatyName_RandomDistribuitionUsed_NatureOfThePotencial(diagonal or off diagonal)_
    //  W=PerturbationValue.txt"
    iofile writter;
    string filename = "";
    
    filename +=  string(typeid(distribution).name());
   
    if(off_diag){
      filename += "_OffDiag";
    }

    filename += "_W=" + to_string(int(W[kk]));

    writter.output_data("outputs/ReturnProb_"+filename+".txt", plt_dt1);
    writter.output_data("outputs/StandDev_"+filename+".txt", plt_dt2);
    cout << endl;
  }

  return 0;
}