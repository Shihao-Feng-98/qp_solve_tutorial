/*
Description: Convex-MPC solver
Email: 13247344844@163.com
Author: Shihao Feng
Update time: 2022-10-24
*/
#ifndef CONVEXMPCSOLVER_H
#define CONVEXMPCSOLVER_H

#include <iostream>
#include <vector>
#include <memory>
using namespace std;

#include "spdlog/spdlog.h"
#include "spdlog/fmt/ostr.h"
#include "spdlog/sinks/rotating_file_sink.h"

#include "OsqpEigen/OsqpEigen.h"
#include <Eigen/Dense>
using namespace Eigen;

#define NUM_DOF 12
#define STATE_DIM 13
#define CONSTRAINT_DIM 20 
#define HORIZON 10 

inline Matrix3d skew(const Vector3d &v)
{
    Matrix3d skew_mat;
    skew_mat << 0., -v(2), v(1),
                v(2), 0., -v(0),
                -v(1), v(0), 0.;
    return skew_mat;
};

class ConvexMpcSolver
{
public:
    ConvexMpcSolver(double dt_mpc, VectorXd Q, double R);
    ~ConvexMpcSolver();
    VectorXd solve(VectorXd body_state, 
                    VectorXd body_state_d, 
                    vector<Vector3d> p_BF_list, 
                    vector<bool> contact_schedule);

private:
    void reset();
    void _calc_A_mat(double yaw_WB);
    // void _calc_contact_scheduled();
    void _calc_B_mat(Vector3d euler_WB, vector<Vector3d> p_BF_list, vector<bool> contact_schedule);
    void _calc_horizon_state(VectorXd body_state, VectorXd body_state_d);
    void _calc_qp_mat(vector<bool> contact_schedule);

    double _dt_mpc; // dt of mpc, may different from the dt of motor control

    double _m; // mass
    DiagonalMatrix<double,3> _I_BB; // inertia matrix w.r.t body frame
    double _mu; // friction coeffient
    double _Fz_min, _Fz_max; // force limit
    Matrix<double, 5, 3> _linear_constraint_leg;

    // QP problem matrices and vectors
    DiagonalMatrix<double, STATE_DIM * HORIZON> _Q_mat; // 130x130
    DiagonalMatrix<double, NUM_DOF * HORIZON> _R_mat; // 120x120

    SparseMatrix<double> _hessian_mat, _linear_constraint_mat; 
    VectorXd _gradient_mat, _lb_mat, _ub_mat; 
    // Matrix<double, NUM_DOF * HORIZON, 1> _gradient_mat; // 120x1
    // Matrix<double, CONSTRAINT_DIM * HORIZON, 1> _lb_mat, _ub_mat; // 200x1
    
    // mpc problem matrices and vectors
    Matrix<double, STATE_DIM, 1> _cur_state; // 13x1
    Matrix<double, STATE_DIM * HORIZON, 1> _mpc_state_d; // 130x1
    vector<bool> _contact_state; 

    Matrix<double, STATE_DIM, STATE_DIM> _Ac_mat, _Ad_mat; // 13x13
    Matrix<double, STATE_DIM, NUM_DOF> _Bc_mat, _Bd_mat; // 13x12
    Matrix<double, STATE_DIM * HORIZON, STATE_DIM> _Aqp_mat; // 130x13
    Matrix<double, STATE_DIM * HORIZON, NUM_DOF * HORIZON> _Bqp_mat; // 130x120

    // // for debug
    // shared_ptr<spdlog::logger> logger;
};


#endif
