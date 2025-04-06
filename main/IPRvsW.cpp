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
    vector<double> L = {12, 32, 64, 124}; // # of lattice points
    const int N = 30; // # of simulations per parameter for the averaging of the simulation values
    vector<double> T_pars = {-5, 20, 0.5, 1}; // Time specifications of the simulation {tmin, tmax, step, log scale(0:false, 1:true)}
    vector<double> W_pars = {0.2, 4, 0.1}; // Perturbation specifications of the simulation {tmin, tmax, step, log scale(0:false, 1:true)}
    vector<double> W; // perturbance amplitude
    for (double uu = W_pars[0]; uu < W_pars[1]; uu += W_pars[2])
    {
        W.push_back(uu);
    }
    bool off_diag = false;

    // Distribution function initialization
    default_random_engine generator;

    // IPR - Eigen values
    vector<vector<double>> plt_dt1; // data for output in .txt file: {tpars, Lvals, data}
    vector<double> data1, aveg_data1;
    plt_dt1.push_back(W_pars);
    plt_dt1.push_back(L);

    // Parallelizing the loop over W
    for (int kk = 0; kk < W.size(); kk++)
    {
        //cout << "W=" << W[kk] << ":" << endl; 
        printProgressBar(kk, W.size());

        // Distribution function:
        uniform_real_distribution<double> distribution(0, 2 * M_PI);

        // Vector to store results for all L values in this W iteration
        vector<double> avg_data_for_W;

        #pragma omp parallel for
        for (int ii = 0; ii < L.size(); ii++)
        {
            //cout << "-L=" << L[ii] << endl;

            // Parameters of the simulation
            vector<int> pars = {int(L[ii]), int(W[kk])};

            // Initialization of the Hamiltonian matrix
            MatrixXcd ham(pars[0], pars[0]);

            vector<double> ipr_per_L;

            // Perform simulations for N realizations
            #pragma omp parallel for
            for (int ll = 0; ll < N; ll++)
            {
                // Hamiltonian creation with a diagonal potential perturbation
                ham_PeridPot_creator(ham, W[kk], distribution(generator));

                // Eigen-solver (use Eigen's ComplexEigenSolver)
                ComplexEigenSolver<MatrixXcd> ces(ham);

                // Extract eigenvectors and compute the IPR (Inverse Participation Rati
                
                MatrixXcd eigen_vcts = ces.eigenvectors();
                VectorXcd aux;
                vector<double> aux2;
                
                eigen_vcts = eigen_vcts.array().pow(4);
                aux = eigen_vcts.rowwise().sum();
                    
                ipr_per_L.push_back(aux.sum().real()/L[ii]);

            }

            // Average over N realizations
            avg_data_for_W.push_back(mean_vector(ipr_per_L));

        }

        // Store results for the current W value
        #pragma omp critical
        {
            plt_dt1.push_back(avg_data_for_W);
        }

    }

    // Output the results to a file
    iofile writter;
    string filename = "Wvs.IPR_QuasiPeriodicPert_withAlpha_";
    
    if (off_diag)
    {
        filename += "_OffDiag";
    }

    writter.output_data("outputs/" + filename + ".txt", plt_dt1);

    return 0;
}
