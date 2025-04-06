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
  cout << "Insert the system's lattice size you wish to analyse:\n";
  vector<double> L = readVectorFromInput(); // # of lattice points
  cout << "Insert the system's perturbation amplitude you wish to analyse:\n";
  vector<double> W = readVectorFromInput(); // perturbation amplitude
  cout << "Random potencial field in each lattice (false) or random hopping potential (true):\n";
  bool off_diag;
  cin >> off_diag; // Parameter which controls if the hamiltonian has an random potential
                   // field at each lattice point (false) or if the hopping potential is
                   // randomized (true)

  const int N = 30; // # of simulations per parameter for the averaging of the simulation values
  vector<double> T_pars = {-5, 20, 0.5, 1}; // Time specifications of the simulation 
                                            // {tmin, tmax, step, log scale(0:false, 1:trues)}
  bool PBC = true;

  // Distribution function initialization
  default_random_engine generator;

  // loops through the desired amplitude perturbation values
  #pragma omp parallel for
  for(int kk = 0; kk<W.size(); kk++)
  {    
    // Define thread-local copies
    vector<double> data;
    vector<vector<double>> all_vals, final_vct;
    VectorXd vaux;

    // Distribution function:
    uniform_real_distribution<double> distribution(0, 2*M_PI);

    // Parameters of the simulation
    vector<int> pars = {int(L[0]), int(W[kk])};
    vector<double> prt = {W[kk]};

    final_vct.push_back(L);
    final_vct.push_back(prt);

    vaux.resize(pars[0]);

    // Innitialization of the Matrix object which will be the hamiltonian
    MatrixXcd ham(pars[0], pars[0]);

    // loops throught the beta values
    for (long double ii=0; ii< 0.3200; ii+=0.0002)
    {
      double beta = 2 * M_PI * ii;

      for (int ll=0; ll<N; ll++)
      {
        printProgressBar(ll, N);

        // Hamiltonian creation 
        ham_PeridPot_creator(ham, W[kk], beta, distribution(generator), off_diag, PBC);

        // Eigen-vector and -values solver
        ComplexEigenSolver<MatrixXcd> ces(ham);
        
        // Storage of the eigenvalues
        vaux = ces.eigenvalues().real();
        eigenvct_to_vector(vaux, data);
        all_vals.push_back(data);

        data.clear();        
      }

      // Averaging of the eigenvalues for the same beta value
      final_vct.push_back(mean_vector(all_vals));
      all_vals.clear();
    }

    // Creating .txt file where the simulation's information will be stored to be plot later
    // with the main.py. The text file name has the following structure:
    // "EnergySpectrum_QuasiPeriodicPert_withAlpha_NatureOfThePotencial(diagonal or off diagonal)_
    //  IfPBCorNot_W=PerturbationValue_L=SystemSizeValue.txt"
    iofile writter;
    string filename = "";
    filename +=  "QuasiPeriodicPert_withAlpha_";

    if(off_diag){
      filename += "_OffDiag";
    }
        
    if(PBC){
      filename += "_PBC";
    }

    writter.output_data("outputs/EnergySpectrum_"+filename+"_W=" + to_string(int(W[kk]))+"_L=" + to_string(int(L[0]))+".txt", final_vct);

    final_vct.clear();
  }

  cout<< endl;

  return 0;
}