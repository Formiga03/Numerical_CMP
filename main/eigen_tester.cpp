#include <iostream>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
#include <vector>
#include <complex>
using namespace std;
using namespace Eigen;
using Complex = std::complex<float>;

void exp_vect_vals(VectorXcd& vect)
{
  for(int ii = 0; ii < vect.size(); ++ii)
  {
    cout << exp(vect(ii)) <<endl;
    vect(ii) = exp(vect(ii));
  }
}

void abs2_vect_vals(VectorXcd& vect)
{
  for(int ii = 0; ii < vect.size(); ++ii)
  {
    vect(ii) = abs(vect(ii))*abs(vect(ii));
  }
}



int main()
{

  vector<double> vct = {1,2};
  for(int ii=0; ii<2; ii++){
    vector<double> vect1;
    for (int jj=0; jj<10; jj++)
    {
      vect1.push_back(vct[ii]*jj);
    }

    cout << vect1.size() << endl;

    cout<<"-------------"<< endl;
  }

  
}


/*

  for (double xx=-5; xx<20; xx=xx+0.5)
  {
    cout << pow(2,xx) << " ";
  }

  cout << endl;



  const double pi = std::acos(-1.0);
  MatrixXcf m1(3,3);

  MatrixXcf m2(3,3);
  m2 = MatrixXcf::Random(3,3);

  MatrixXcf m3(3,3);
  m3(0,1) = Complex(pi/4, 0);
  m3(1,0) = Complex(-pi/4, 0);

  


  ComplexEigenSolver<MatrixXcf> ces(m3);


  //cout << ces.eigenvalues() << endl <<endl;
  //cout << ces.eigenvectors().exp()<< endl;

  VectorXcd vect = VectorXcd::Zero(3);
  VectorXcd vect2(vect);

  vect2 << Complex(1,2),
            Complex(1,0),
            Complex(0,0);

  //ComplexScalar<double> z1 = (1,2);
  //vect2(1)= (1,0);

  //cout << vect2*complex<double>(0,-1) << endl;

  for(int i = 0; i < vect2.size(); ++i)
  {
      cout << exp(vect2(i)) << endl;
  }

  complex<double> z1(-3.49854778,0);

  complex<double> res = exp(z1*complex<double>(0,-1));

  cout << z1 << " | " << res << endl;

*/