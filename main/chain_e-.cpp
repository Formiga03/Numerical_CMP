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
  bool even=true;
  const int N = 30; // # of simulations per parameter for the averaging of the simulation 
                    // values
  vector<double> T_pars = {-5, 20, 0.5, 1}; // Time specifications of the simulation
                                            // {tmin, tmax, step, log scale(0:false, 1:trues)}
    
  
  // Distribution function initialization
  default_random_engine generator;

  for(int kk = 0; kk<W.size(); kk++)
  {
    cout << "W=" << W[kk] << ":" << endl; 

    // Distribution function:
    normal_distribution<double> distribution(2.0,W[kk]);

    // distribuition function and perturbation info
    vector<double> pert = {W[kk]};

    // imbalance
    vector<vector<double>> plt_dt1; // data for output in .txt file: {tpars, Lvals, data}
    plt_dt1.push_back(pert);
    plt_dt1.push_back(T_pars);
    plt_dt1.push_back(L);

    // entangement entropy
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
      MatrixXcd rho_0(pars[0], pars[0]);
      MatrixXcd phi_0(pars[0], pars[0]/2);
      phi_0.setZero();

      // Initialization of the density matrix where we will store the evolution in time of the value theat characterize every simulated system
      vector<vector<double>> sim_data1, sim_data2;

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

        
        VectorXcd aux(pars[0]/2);
        if (even){
          for (int jjj=0; jjj<pars[0]; jjj++)
          {
            aux.setZero();
            if (jjj%2 == 1){
              aux((jjj-1)/2) = 1;
              phi_0.row(jjj) = aux;
            }
          }
        } if (not even){
          for (int jjj=0; jjj<pars[0]; jjj++)
          {
            aux.setZero();
            if (jjj%2 == 0){
              aux(jjj/2) = 1;
              phi_0.row(jjj) = aux;
            }
          }
        }

        MatrixXcd rho(rho_0);
        MatrixXcd phi(pars[0], pars[0]);
        phi.setZero();
    
        VectorXcd eigen_vls(ces.eigenvalues());
        VectorXcd aux_vct1(pars[0]);
        VectorXcd aux_vct2(pars[0]);
        VectorXcd aux_vct3(pars[0]);
        VectorXcd aux_vct4(pars[0]);
        VectorXcd rho_diag(pars[0]);
        MatrixXcd eigen_vcts(ces.eigenvectors());

        MatrixXcd aux_mat1(pars[0], pars[0]);
        MatrixXcd aux_mat2(pars[0], pars[0]);
        MatrixXcd aux_mat3(pars[0], pars[0]);
        MatrixXcd aux_mat4(pars[0], pars[0]);

        MatrixXcd aux_mat5(pars[0]/2, pars[0]/2);
        
        vector<double> data1, data2;

        // Initial value of time
        double tt;

        for (double xx=T_pars[0]; xx<T_pars[1]; xx+=T_pars[2])
        {
          tt = pow(2, xx); 
          
          aux_vct1.setZero();
          aux_vct2.setZero();
          aux_vct3.setZero();
          aux_vct4.setZero();
          rho_diag.setZero();

          aux_mat1.setZero();
          aux_mat2.setZero();
          aux_mat3.setZero();
          aux_mat4.setZero();
          aux_mat5.setZero();

          for(int jjj = 0; jjj < pars[0]; ++jjj)
          {
            aux_vct1(jjj) = eigen_vls(jjj) * complex<double>(0,-1) * complex<double>(tt,0);
            //aux_vct2(jjj) = eigen_vls(jjj) * complex<double>(0,1) * complex<double>(tt,0);
          }

          exp_vect_vals(aux_vct1);
          //exp_vect_vals(aux_vct2);

          aux_mat1 = aux_vct1.asDiagonal(); // diagonalized time evolution exponencial matrix
          //aux_mat2 = aux_vct2.asDiagonal(); // hermitian diagonalized time evolution exponencial matrix

          aux_mat3 = eigen_vcts * aux_mat1 * eigen_vcts.conjugate().transpose(); // time evolution mat. 
          //aux_mat4 = eigen_vcts * aux_mat2 * eigen_vcts.conjugate().transpose(); // hermitian time evolution mat. 

          rho = aux_mat3 * rho_0 * aux_mat3.conjugate().transpose(); 
          rho_diag = rho.diagonal();

          complex<double> even_ind = 0; 
          complex<double> odd_ind = 0; 

          for (int ee = 0; ee < rho_diag.size(); ee += 2) 
          {
            even_ind += rho_diag(ee+1);
            odd_ind += rho_diag(ee);
          }

          if (even) {
            data1.push_back(2*even_ind.real()/pars[0] - 2*odd_ind.real()/pars[0]);
          } if (not even) {
            data1.push_back(2*odd_ind.real()/pars[0] - 2*even_ind.real()/pars[0]);
          }

          // Entanglement Entropy Calc.
          
          phi = aux_mat3 * phi_0 * phi_0.conjugate().transpose() * aux_mat3.conjugate().transpose();
          aux_mat5 = phi.block(0, 0, pars[0] / 2, pars[0] / 2);

          ComplexEigenSolver<MatrixXcd> ces1(aux_mat5);

          double threshold = pow(10, -20);

          // Compare using absolute values
          aux_vct3 = (ces1.eigenvalues().array().abs() >= threshold).select(ces1.eigenvalues(), 0);
          aux_vct4 = VectorXcd::Ones(pars[0] / 2) - ces1.eigenvalues();
          aux_vct4 = (aux_vct4.array().abs() >= threshold).select(aux_vct4, 0);

          complex<double> S = -aux_vct3.dot(log_vect_vals(aux_vct3)) - aux_vct4.dot(log_vect_vals(aux_vct4));
          data2.push_back(S.real());

        }

        sim_data1.push_back(data1);
        sim_data2.push_back(data2);
      }
      cout<< endl;

      vector<double> graph_data1, graph_data2;

      mean_vector(sim_data1, graph_data1);
      mean_vector(sim_data2, graph_data2);
      
      plt_dt1.push_back(graph_data1);
      plt_dt2.push_back(graph_data2);

    }
    cout<< endl;
    iofile writter;
    string filename =  string(typeid(distribution).name());
    if(off_diag){
      filename += "_OffDiag";
    }
    filename += "_W=" + to_string(int(W[kk]));

    writter.output_data("outputs/Imbalance_"+filename+".txt", plt_dt1);
    writter.output_data("outputs/EntenglEntrop_"+filename+".txt", plt_dt2);
  }
  cout<< endl;

  return 0;
}