#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <complex>
#include <random>
#include <cmath>


using namespace std;
using namespace Eigen;

void eigenvct_to_vector(VectorXd vct1, vector<double> vct2);
void print_ham_eigen_vectors_values(MatrixXcd ham, ComplexEigenSolver<MatrixXcd> ces, vector<int> pars);
void exp_vect_vals(VectorXcd& vect);
VectorXcd log_vect_vals(VectorXcd& vect);
void abs2_vect_vals(VectorXcd& vect);
void mean_vector(vector<vector<double>>& vct, vector<double>& res);
void hamiltonian_creator(MatrixXcd& ham,default_random_engine& gen, normal_distribution<double>& dist, bool off_diag=false);
void hamiltonian_creator(MatrixXcd& ham, default_random_engine gen, uniform_real_distribution<double> dist, bool off_diag=false); 
void hamiltonian_creator(MatrixXcd& ham, default_random_engine gen, gamma_distribution<double> dist, bool off_diag=false);
void ham_PeridPot_creator(MatrixXcd& ham, double W, double aplha, bool off_diag=false);
void printProgressBar(int progress, int total, int barWidth = 50);
void density_creator(MatrixXcd& dens, bool even=true);


