#include <iostream>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
using namespace std;
using namespace Eigen;

void print_mat(MatrixXf mat){
    int rows = mat.rows();
    int cols = mat.cols();

    for (int ii=0; ii<rows; ii++)
    {
        for (int jj=0; jj<cols; jj++)
        {
            cout << mat(ii, jj) << " ";
        }
        cout << "\n";
    }
    cout << endl;

}


int main(){
    const double pi = std::acos(-1.0);
    MatrixXf m1(3,3);

    MatrixXf m2(3,3);
    m2 = MatrixXf::Random(3,3);

    MatrixXf m3(3,3);
    m3 << 0,    -pi/4, 0,
       pi/4, 0,     0,
       0,    0,     0;

    cout << m3 << endl <<endl;
    cout << m3.exp() << endl;

    VectorXcd vect=VectorXcd::Zero(3);
    VectorXcd vect2(vect);

    vect2(1)=(1, 0);
    cout << vect2*(0,-1) << endl;
    cout << vect << endl;
    

    //print_mat(m1);
    //print_mat(m2);
    //print_mat(m3);


}
