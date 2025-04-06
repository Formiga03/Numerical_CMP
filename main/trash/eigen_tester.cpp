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
#include <Eigen/Sparse>
#include <Eigen/Eigenvalues>

//#include "../inc/matplotlibcpp.h"
#include "../inc/functions.h"
#include "../inc/iofile.h"

using namespace std;
using namespace Eigen;
//namespace plt = matplotlibcpp;

void print3DVector(const vector<vector<vector<double>>>& vec) {
    for (size_t i = 0; i < vec.size(); ++i) {
        cout << "Layer " << i << ":\n";
        for (size_t j = 0; j < vec[i].size(); ++j) {
            for (size_t k = 0; k < vec[i][j].size(); ++k) {
                cout << vec[i][j][k] << " ";
            }
            cout << endl;
        }
        cout << endl;
    }
}
void printVector(const std::vector<double>& vec) {
    for (const auto& val : vec) {
        std::cout << val << " ";
    }
    std::cout << std::endl;
}

int main() 
{

    std::cout << "Enter values for vector 1: ";
    std::vector<double> vec1 = readVectorFromInput();

    std::cout << "Vector 1: ";
    for (auto& n : vec1) std::cout << n << " ";
    



    /*
    vector<int> L = {4};

    SparseMatrix<double> mat(3, 3);
    mat.insert(0, 0) = 4;
    mat.insert(1, 1) = 4;
    mat.insert(2, 2) = 4;
    mat.insert(0, 1) = -2;
    mat.insert(1, 0) = -2;
    mat.insert(1, 2) = -2;
    mat.insert(2, 1) = -2;
    mat.insert(0, 2) = 1;
    mat.insert(2, 0) = 1;


    for (auto& ee: L) {
       
        MatrixXcd ham(ee, ee);  // Each thread gets its own matrix

        //#pragma omp parallel for
        for (int ii = 0; ii < 1; ii++) {
            //printProgressBar(ii,30);
            cout << endl;
            ham_PeridPot_creator(ham, 4, 0.5, false, true);
            
            ComplexEigenSolver<MatrixXcd> solver(ham);
            VectorXcd vct = solver.eigenvalues();

            cout << vct << endl;
            cout << "_____________________" << endl;
            cout << vct.array().pow(2) << endl;
            cout << "_____________________" << endl;
            cout << vct.array().exp() << endl;
            cout << "_____________________" << endl;
            cout << vct.cwiseAbs().array().log() << endl;
            cout << "_____________________" << endl;
        }
        cout << endl;

    }



/*

    Eigen::VectorXd vctaux(3);
    mat << 1, 2, 3,
           1, 2, 3,
           0, 1, 2;
    mat.normalize();

    vector<double> aux;
    double threshold = 2.0;
    
    // Apply the thresholding operation
    //mat = (mat.array() >= pow(10,-20)).select(mat, 0);
    vctaux = mat.colwise().norm().array().pow(2).real();
    std::cout << mat << std::endl;
    std::cout << vctaux  << std::endl;
    eigenvct_to_vector(vctaux, aux);
    printVector(aux);

    vector<vector<vector<double>>> vct = {{{1}, {3}}, {{2}, {4}}};
    print3DVector(vct);
    W_sort(vct);
    print3DVector(vct);

    int n = 4;
    MatrixXcd phi_0(n,n/2);
    VectorXcd aux(n/2);
    for (int jjj=0; jjj<n; jjj++)
    {
        aux.setZero();
        if (jjj%2 == 1){
            aux((jjj-1)/2) = 1;
            cout << jjj << " | " << (jjj+1)/2 << " | " << aux.transpose() << endl ;
            phi_0.row(jjj) = aux;
        }
    }
    cout << phi_0* phi_0.transpose() << endl;

    MatrixXcd psi_0(n,n);
    density_creator(psi_0);

    default_random_engine generator;
    normal_distribution<double> distribution(2.0,2);
    MatrixXcd ham(n,n);
    //ham_PeridPot_creator(ham);
    cout << ham << endl;

    for (int i =0; i<1; i++)
    {
        ham.setZero();
        hamiltonian_creator(ham, generator, distribution);
        cout << ham << endl << endl;
        cout << 2*ham/n << endl << endl;
    }
    */




}
