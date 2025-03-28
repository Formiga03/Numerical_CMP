#include "functions.h"

void eigenvct_to_vector(VectorXd vct1, vector<double> vct2)
{
  for (int i = 0; i < vct1.size(); ++i) 
  {
    vct2.push_back(vct1(i));    
  }
}

void print_ham_eigen_vectors_values(MatrixXcd ham, ComplexEigenSolver<MatrixXcd> ces, vector<int> pars)
{
  cout << "Hamiltonian: (L=" << pars[0] << ", W=" << pars[1] << ", Iter=" << pars[2] << ")" << endl;
  cout << ham << endl;
  cout << "_______________________________________________" << endl;
  cout << "Eigenvectors:" << endl;
  cout << ces.eigenvectors() << endl;
  cout << "_______________________________________________" << endl;
  cout << "Eigenvalues:" << endl;
  cout << ces.eigenvalues() << endl;
  cout << "_______________________________________________\n" << endl;

}

void exp_vect_vals(VectorXcd& vect)
{
  for(int ii = 0; ii < vect.size(); ++ii)
  {
    vect(ii) = exp(vect(ii));
  }
}

VectorXcd log_vect_vals(VectorXcd& vect)
{
  VectorXcd res(vect.size());
  for(int ii = 0; ii < vect.size(); ++ii)
  {
    res(ii) = log(vect(ii));
  }

  return res;
}

void abs2_vect_vals(VectorXcd& vect)
{
  MatrixXcd aux = vect.conjugate().asDiagonal();
  vect = aux * vect;
}

void mean_vector(vector<vector<double>>& vct, vector<double>& res)
{
  double aux_mean;

  for(int ii=0; ii<vct[0].size(); ii++)
  {
    aux_mean = 0;
    for(int jj=0; jj<vct.size(); jj++)
    {
      aux_mean += vct[jj][ii];
    }

    res.push_back(aux_mean/vct.size());
  }
}

void hamiltonian_creator(MatrixXcd& ham,default_random_engine& gen, normal_distribution<double>& dist, bool off_diag) 
{
  if (off_diag) {
    double aux_num;
    for (int jj=0; jj<ham.rows()-1; jj++)
    {
      aux_num = dist(gen);
      ham(jj, jj+1) = aux_num; 
      ham(jj+1, jj) = aux_num; 
    }
  } else {
    for (int jj=0; jj<ham.rows()-1; jj++)
    {
      ham(jj, jj) = dist(gen);
      ham(jj, jj+1) = 1; 
      ham(jj+1, jj) = 1; 
    }
    ham(ham.rows()-1,ham.rows()-1) = dist(gen);
  }
}

void hamiltonian_creator(MatrixXcd& ham, default_random_engine& gen, uniform_real_distribution<double>& dist, bool off_diag) 
{
  if (off_diag) {
    for (int jj=0; jj<ham.rows()-1; jj++)
    {
      ham(jj, jj+1) = dist(gen); 
      ham(jj+1, jj) = dist(gen); 
    }
  } else {
    for (int jj=0; jj<ham.rows()-1; jj++)
    {
      ham(jj, jj) = dist(gen);
      ham(jj, jj+1) = 1; 
      ham(jj+1, jj) = 1; 
    }
    ham(ham.rows()-1,ham.rows()-1) = dist(gen);
  }
}

void hamiltonian_creator(MatrixXcd& ham, default_random_engine& gen, gamma_distribution<double>& dist, bool off_diag) 
{
  if (off_diag) {
    for (int jj=0; jj<ham.rows()-1; jj++)
    {
      ham(jj, jj+1) = dist(gen); 
      ham(jj+1, jj) = dist(gen); 
    }
  } else {
    for (int jj=0; jj<ham.rows()-1; jj++)
    {
      ham(jj, jj) = dist(gen);
      ham(jj, jj+1) = 1; 
      ham(jj+1, jj) = 1; 
    }
    ham(ham.rows()-1,ham.rows()-1) = dist(gen);
  }
}

void ham_PeridPot_creator(MatrixXcd& ham, double W, double alpha, bool off_diag)
{
  ham.setZero();
  const double sgm = (1+sqrt(5))/2;
  auto u = [sgm, W, alpha](double ii) -> double {
    double res = W*cos(2*M_PI*sgm*ii + alpha);
    return res;
  };

  for (int iii=0; iii<ham.rows(); iii++)
  {
    ham(iii, iii) = u(iii);
  }
  for (int jj=0; jj<ham.rows()-1; jj++)
  {
    ham(jj, jj+1) = 1; 
    ham(jj+1, jj) = 1;
  }
}

void print_vals(vector<double> data1, vector<double> data2)
{
  for (int oo=0; oo<data1.size(); oo++)
    {
      cout << oo << ": "<<data1[oo] << " | " << data2[oo] << endl;
    }
}

void printProgressBar(int progress, int total, int barWidth) 
{
  float percent = (float)progress / total;
  int filled = percent * barWidth;

  std::cout << "\r[";
  for (int i = 0; i <= barWidth; i++) {
      if (i < filled)
          std::cout << "=";
      else if (i == filled)
          std::cout << ">";
      else
          std::cout << " ";
  }
  std::cout << "] " << int(percent * 100.0) << "%";
  std::cout.flush();
}

void density_creator(MatrixXcd& dens, bool even)
{
  int dim =  dens.rows();

  dens.setZero();
  Eigen::VectorXd diag_elements(dim);

  for (int hh = 0; hh < dim; ++hh) 
  {
    diag_elements(hh) = ((hh) % 2 == 1) ? 1 : 0; 
  }

  dens = diag_elements.asDiagonal();
}
