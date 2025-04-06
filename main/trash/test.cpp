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
    vector<double> L = {986}; // # of lattice points
    vector<double> W = {1, 2}; // perturbance amplitude
    const int N = 30; // # of simulations per parameter
    bool PBC = true;
    bool off_diag = false;

    vector<vector<double>> final_vct;
    final_vct.reserve(W.size() * (1 + 1 / 0.01) * N); // Preallocate memory

    // Distribution function initialization
    default_random_engine generator;
    uniform_real_distribution<double> distribution(0, 2 * M_PI);

    #pragma omp parallel for schedule(dynamic)
    for (size_t kk = 0; kk < W.size(); kk++)
    {
        vector<vector<double>> all_vals;
        all_vals.reserve(N);

        cout << "W=" << W[kk] << ":" << endl; 
        vector<int> pars = {int(L[0]), int(W[kk])};
        vector<double> prt = {W[kk]};

        vector<double> data;
        data.reserve(L[0]);
        
        final_vct.push_back(L);
        final_vct.push_back(prt);

        MatrixXcd ham(pars[0], pars[0]);
        
        for (double ii = 0; ii <= 1; ii += 0.01)
        {
            cout << ii << endl;
            double beta = ii * (1 + sqrt(5)) / 2;

            #pragma omp parallel for
            for (int ll = 0; ll < N; ll++)
            {
                MatrixXcd ham_local = ham; // Thread-safe local copy
                ham_PeridPot_creator(ham_local, W[kk], beta, distribution(generator), off_diag, PBC);
                
                ComplexEigenSolver<MatrixXcd> ces(ham_local);
                
                vector<double> eigenvalues(pars[0]);
                for (int i = 0; i < pars[0]; i++) {
                    eigenvalues[i] = ces.eigenvalues()[i].real();
                }
                
                #pragma omp critical
                all_vals.push_back(eigenvalues);
            }
            
            #pragma omp critical
            final_vct.push_back(mean_vector(all_vals));
            all_vals.clear();
        }
    }

    // Writing output file after computation
    iofile writter;
    string filename = "QuasiPeriodicPert_withAlpha";
    if (PBC) filename += "_PBC";
    if (off_diag) filename += "_OffDiag";

    for (size_t kk = 0; kk < W.size(); kk++) {
        writter.output_data("outputs/EnergySpectrum_" + filename + "_W=" + to_string(int(W[kk])) + "_L=" + to_string(int(L[0])) + ".txt", final_vct);
    }
    
    cout << "Computation complete." << endl;
    return 0;
}