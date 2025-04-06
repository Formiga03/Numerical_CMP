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
#include <omp.h>

#include "../inc/functions.h"
#include "../inc/iofile.h"

using namespace std;
using namespace Eigen;

int main() 
{
    // Parameters
    vector<int> L = {200}; // Lattice size
    vector<double> W = {1.0, 2.0}; // Perturbation amplitude
    const int N = 30; // Number of simulations per parameter for averaging
    const bool PBC = true;
    const bool off_diag = false;

    vector<double> data;
    vector<vector<double>> final_vct;

    // Random number generator setup
    random_device rd;

    for (size_t kk = 0; kk < W.size(); kk++)
    {
        cout << "W = " << W[kk] << ":" << endl;

        vector<int> pars = {L[0], static_cast<int>(W[kk])};
        vector<double> prt = {W[kk]};
        
        final_vct.push_back({static_cast<double>(L[0])});
        final_vct.push_back(prt);

        VectorXd vaux(pars[0]);
        MatrixXcd ham(pars[0], pars[0]);

        // Parallelize only the inner loops
        #pragma omp parallel 
        {
            default_random_engine generator(rd()); // Thread-local generator
            uniform_real_distribution<double> distribution(0, 2 * M_PI);

            #pragma omp for schedule(dynamic) collapse(2)
            for (int ii = 0; ii <= pars[0]; ii++)
            {
                //printProgressBar(ii, pars[0]);
                double beta = static_cast<double>(ii) / pars[0];

                vector<vector<double>> all_vals; // Thread-local storage

                for (int ll = 0; ll < N; ll++)
                {
                    // Create Hamiltonian with perturbation
                    ham_PeridPot_creator(ham, W[kk], beta, distribution(generator), off_diag, PBC);

                    // Solve for eigenvalues
                    ComplexEigenSolver<MatrixXcd> ces(ham);
                    vaux = ces.eigenvalues().real();

                    // Store results
                    vector<double> local_data;
                    eigenvct_to_vector(vaux, local_data);
                    all_vals.push_back(local_data);
                }

                // Merge results in a thread-safe manner
                #pragma omp critical
                {
                    final_vct.push_back(mean_vector(all_vals));
                }
            }
        }

        // File Writing
        iofile writer;
        string filename = "outputs/EnergySpectrum_QuasiPeriodicPert_withAlpha_";

        if (PBC) filename += "_PBC";
        if (off_diag) filename += "_OffDiag";

        filename += "_W=" + to_string(static_cast<int>(W[kk])) + "_L=" + to_string(L[0]) + ".txt";

        writer.output_data(filename, final_vct);

        final_vct.clear();
    }

    cout << endl;
    return 0;
}
