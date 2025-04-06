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
  cout << "Periodic boundary conditions (true/false):\n";
  bool PBC;
  cin >> PBC;

  const int N = 30; // # of simulations per parameter for the averaging of the simulation values
  vector<double> T_pars = {-5, 20, 0.5, 1}; // Time specifications of the simulation {tmin, tmax, step, log scale(0:false, 1:trues)}
  

  // Distribution function initialization
  default_random_engine generator;

  // data for output in .txt file: {Lval, tpars, Wvals, data}
  vector<vector<double>> plt_dt4;
  vector<vector<vector<double>>> GrdSt_dt;

  // Initialization of the vectors where we will store the evolution in time of the value theat characterize every simulated system
  vector<vector<double>> sim_data4;
  vector<double> data4;
  vector<double> graph_data4;

  // Parameters of the simulation
  vector<int> pars;


  for(int kk = 0; kk<W.size(); kk++)
  {
    cout << "W=" << W[kk] << ":" << endl; 

    // Distribution function:
    uniform_real_distribution<double> distribution(0, 2*M_PI);

    for(int ii = 0; ii<L.size(); ii++)
    {
      cout << "-L=" << L[ii] << endl;
    
      pars = {int(L[ii]), int(W[kk])};

      // Innitialization of the Matrix object which will be the hamiltonian
      MatrixXcd ham(pars[0], pars[0]);

      #pragma omp parallel for
      for (int ll=0; ll<10; ll++)
      {
        printProgressBar(ll, N);

        // Hamiltonian creation with a diagonal potential perturbation
        ham_PeridPot_creator(ham, W[kk], distribution(generator), off_diag, PBC);
        
        // Eigen-vector and -values solver
        ComplexEigenSolver<MatrixXcd> ces(ham);
        VectorXd eigen_vls(ces.eigenvalues().real());
        MatrixXcd eigen_vcts(ces.eigenvectors());

        // Ground State Calc.:
        VectorXd aux_vct(pars[0]);
        double maxValue = eigen_vls.minCoeff();
        int maxIndex;
        eigen_vls.minCoeff(&maxIndex);
        aux_vct = eigen_vcts.col(maxIndex).cwiseAbs();     
        eigenvct_to_vector(aux_vct, data4); 
        sim_data4.push_back(data4); 

        data4.clear();
        cout << endl;

      }
      cout << endl;

      // Averaging for the number of hamiltonians with different random phase
      mean_vector(sim_data4, graph_data4); 

      plt_dt4.push_back(graph_data4);
      
      sim_data4.clear();
      graph_data4.clear();

    }
    cout<< endl;
 
    GrdSt_dt.push_back(plt_dt4);
    
    plt_dt4.clear();

  }

  iofile writter;
  string filename = "";
  filename +=  "QuasiPeriodicPert_withAlpha_W={";

  for (auto uu: W)
  {
    filename += to_string(int(uu)) +" ";
  }

  filename += "\r}_";

  vector<vector<double>> aux_vect;
  vector<double> aux_vect1;

  W_sort(GrdSt_dt);

  if (PBC){
      filename += "PBC";
  }

  for (int kkk=0; kkk<L.size(); kkk++)
  {
    aux_vect1.push_back(L[kkk]);

    if (PBC){
      aux_vect1.push_back(1);
    } else {
      aux_vect1.push_back(0);
    }

    aux_vect.push_back(aux_vect1);
    aux_vect.push_back(W);

    for (auto ee: GrdSt_dt[kkk])
    {
      aux_vect.push_back(ee);
    }

    writter.output_data("outputs/GroundState_"+filename+"_L=" + to_string(int(L[kkk]))+".txt", aux_vect);

    aux_vect.clear();
    aux_vect1.clear();
  }
 
  cout<< endl;

  return 0;
}