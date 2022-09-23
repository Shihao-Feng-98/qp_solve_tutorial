/*
Description: BalanceQP controller solver
Email: 13247344844@163.com
Author: Shihao Feng
Update time: 2022-09-13
*/

#ifndef BALANCE_QP_SOLVER
#define BALANCE_QP_SOLVER

#include <iostream>
using namespace std;
// osqp-eigen
#include "OsqpEigen/OsqpEigen.h"
#include <Eigen/Dense>
using namespace Eigen;


class BalanceQpControllerSolver
{
public:
    BalanceQpControllerSolver(double m, Vector3d I, 
                            double mu, double Fz_min, double Fz_max, 
                            VectorXd Q, double R);
    ~BalanceQpControllerSolver();
    VectorXd solve(bool const *scheduled_contact, 
                    const MatrixXd &A_mat, 
                    const VectorXd &b_mat);

private:
    double _m; // 机身重量
    DiagonalMatrix<double,3> _I_BB; // 机身惯量
    
    double _mu; // 摩擦系数
    double _Fz_min, _Fz_max; // Z方向的力约束
    
    // QP problem matrices and vectors
    DiagonalMatrix<double, 6> _Q_mat;
    MatrixXd _R_mat;
    SparseMatrix<double> _hessian_mat;
    VectorXd _gradient_mat;
    SparseMatrix<double> _linear_mat;
    VectorXd _lb_mat, _ub_mat;
};


#endif
