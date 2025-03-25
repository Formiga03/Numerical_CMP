#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <complex>
#include <random>


using namespace std;
using namespace Eigen;

void eigenvct_to_vector(VectorXd vct1, vector<double> vct2);
void print_ham_eigen_vectors_values(MatrixXcd ham, ComplexEigenSolver<MatrixXcd> ces, vector<int> pars);
void exp_vect_vals(VectorXcd& vect);
void abs2_vect_vals(VectorXcd& vect);
void mean_vector(vector<vector<double>>& vct, vector<double>& res);
void hamiltonian_creator(MatrixXcd& ham,default_random_engine& gen, normal_distribution<double>& dist, bool off_diag=false);
void hamiltonian_creator(MatrixXcd& ham, default_random_engine& gen, uniform_real_distribution<double>& dist, bool off_diag=false); 
void hamiltonian_creator(MatrixXcd& ham, default_random_engine& gen, gamma_distribution<double>& dist, bool off_diag=false);
void printProgressBar(int progress, int total, int barWidth = 50);


