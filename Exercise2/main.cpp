#include <iostream>
#include "Eigen/Eigen"
#include <iomanip>

using namespace std;
using namespace Eigen;

//faccio le dichiarazioni

void ErrorEstimate(const MatrixXd& A,
                  const VectorXd& b,
                  const VectorXd& sol,
                  double& relErrPALU,
                  double& relErrQR);

VectorXd PALUSolver(const MatrixXd& A,
                    const VectorXd& b);

VectorXd QRSolver(const MatrixXd& A,
                    const VectorXd& b);

int main()
{
    Vector2d sol (-1.0000e+0, -1.0000e+0);

    Matrix2d A0 {{5.547001962252291e-01, -3.770900990025203e-02},{8.320502943378437e-01,-9.992887623566787e-01}};

    Vector2d b0 = {-5.169911863249772e-01, 1.672384680188350e-01};

    double relErr0PALU, relErr0QR;
    ErrorEstimate(A0, b0,sol,relErr0PALU, relErr0QR);
    cout<< "1st system: \n";
    cout<< scientific<<setprecision(16)<<"-Solutions-\n PALU: ("<<PALUSolver(A0,b0).transpose()<<")'"<<endl<< "QR: ("<<QRSolver(A0,b0).transpose()<<")'"<<endl;
    cout<< scientific<<setprecision(16)<< "-Relative Errors-\n "<< "PALU: " << relErr0PALU<< " QR: "<<relErr0QR<<endl;

    Matrix2d A1 {{5.547001962252291e-01, -5.540607316466765e-01},{8.320502943378437e-01,-8.324762492991313e-01}};

    Vector2d b1 = {-6.394645785530173e-04, 4.259549612877223e-04};

    double relErr1PALU, relErr1QR;
    ErrorEstimate(A1, b1,sol,relErr1PALU, relErr1QR);
    cout<< "2nd system: \n";
    cout<< scientific<<setprecision(16)<<"-Solutions-\n PALU: ("<<PALUSolver(A1,b1).transpose()<<")'"<<endl<< "QR: ("<<QRSolver(A1,b1).transpose()<<")'"<<endl;
    cout<< scientific<<setprecision(16)<< "-Relative Errors-\n "<< "PALU: " << relErr1PALU<< " QR: "<<relErr1QR<<endl;

    Matrix2d A2 {{5.547001962252291e-01, -5.547001955851905e-01},{8.320502943378437e-01, -8.320502947645361e-01}};

    Vector2d b2 = {-6.400391328043042e-10, 4.266924591433963e-10};

    double relErr2PALU, relErr2QR;
    ErrorEstimate(A2, b2,sol,relErr2PALU, relErr2QR);
    cout<< "3rd system: \n";
    cout<< scientific<<setprecision(16)<<"-Solutions-\n PALU: ("<<PALUSolver(A2,b2).transpose()<<")'"<<endl<< "QR: ("<<QRSolver(A2,b2).transpose()<<")'"<<endl;
    cout<< scientific<<setprecision(16)<< "-Relative Errors-\n "<< "PALU: " << relErr2PALU<< " QR: "<<relErr2QR<<endl;

  return 0;
}

//definisco le funzioni
VectorXd PALUSolver(const MatrixXd& A,
                    const VectorXd& b)
{
    VectorXd solutionPALU = A.fullPivLu().solve(b);
    return solutionPALU;
}

VectorXd QRSolver(const MatrixXd& A,
                    const VectorXd& b)
{
    VectorXd solutionQR = A.colPivHouseholderQr().solve(b);
    return solutionQR;
}

void ErrorEstimate(const MatrixXd& A,
                   const VectorXd& b,
                   const VectorXd& sol,
                   double& relErrPALU,
                   double& relErrQR)
//dando gli errori in referenza posso avere un doppio output
{
    relErrPALU = (PALUSolver(A,b) - sol).norm()/sol.norm();
    relErrQR = (QRSolver(A,b) - sol).norm()/sol.norm();

}

