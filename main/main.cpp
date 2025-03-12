#include <iostream>
#include <Eigen/Dense>
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
    MatrixXf m1(3,3);

    MatrixXf m2(3,3);
    m2 = MatrixXf::Random(3,3);

    MatrixXf m3(3,3);
    m3 = m1*m2;

    cout << m3 << endl;
    //print_mat(m1);
    //print_mat(m2);
    //print_mat(m3);


}
