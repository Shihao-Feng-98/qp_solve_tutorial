#include <ConvexMpcSolver.h>

ConvexMpcSolver::ConvexMpcSolver(double dt_mpc, VectorXd Q, double R)
{
    // // file logger
    // logger = spdlog::rotating_logger_mt("logger", "../logs/rotating_log.txt", 1048576 * 5, 1);
    // logger->set_level(spdlog::level::debug);

    _dt_mpc = dt_mpc;

    _m = 12.6;
    _I_BB.diagonal() << 0.133872, 0.450659, 0.469184;
    _mu = 0.4;
    _Fz_min = 2.0;
    _Fz_max = 150.0;
    _linear_constraint_leg << -1., 0., _mu,
                                1., 0., _mu,
                                0., -1., _mu,
                                0., 1., _mu,
                                0., 0., 1.;

    Matrix<double, STATE_DIM * HORIZON, 1> Q_mpc;
    Matrix<double, NUM_DOF * HORIZON, 1> R_mpc;
    for (int i = 0; i < HORIZON; i++)
    {
        Q_mpc.segment(i*STATE_DIM, STATE_DIM) = Q;
    }
    _Q_mat.diagonal() = Q_mpc;
    R_mpc.setOnes();
    _R_mat.diagonal() = R * R_mpc;

    _hessian_mat.resize(NUM_DOF * HORIZON, NUM_DOF * HORIZON); // 120x120
    _linear_constraint_mat.resize(CONSTRAINT_DIM * HORIZON, NUM_DOF * HORIZON); // 200x120
    _gradient_mat.resize(NUM_DOF * HORIZON, 1); // 120x1
    _lb_mat.resize(CONSTRAINT_DIM * HORIZON, 1); // 200x1
    _ub_mat.resize(CONSTRAINT_DIM * HORIZON, 1); // 200x1
}

ConvexMpcSolver::~ConvexMpcSolver()
{
    // // split the log file
    // logger->warn("\n\n\n");
}

void ConvexMpcSolver::reset()
{
    // _hessian_mat.setZero();
    _gradient_mat.setZero();
    _linear_constraint_mat.setZero();
    _lb_mat.setZero();
    _ub_mat.setConstant(OsqpEigen::INFTY);

    _Ac_mat.setZero();
    _Bc_mat.setZero();
    _Ad_mat.setZero();
    _Bd_mat.setZero();
    _Aqp_mat.setZero();
    _Bqp_mat.setZero();
}

void ConvexMpcSolver::_calc_A_mat(double yaw_WB)
{
    // Theoretically, "horizon" A matrices should be computed 
    // which need to know the actual orientation at each dt_mpc (approximated by the deisred orientation)
    // In engineering, just compute A matrix once
    Matrix3d Rz_WB;
    Rz_WB = AngleAxisd(yaw_WB, Vector3d::UnitZ());
    _Ac_mat.block<3,3>(0,6) = Rz_WB.transpose();
    _Ac_mat.block<3,3>(3,9) = Matrix3d::Identity();
    _Ac_mat(11,12) = 1.;
    // discretization
    _Ad_mat = Matrix<double, STATE_DIM, STATE_DIM>::Identity() + _dt_mpc * _Ac_mat;
}

void ConvexMpcSolver::_calc_B_mat(Vector3d euler_WB, vector<Vector3d> p_BF_list, 
                                    vector<bool> contact_schedule)
{
    // Theoretically, "horizon" B matrices should be computed 
    // which need to know the actual orientation (approximated by the deisred orientation),
    // the contact schedule (approximated by the desired contact schedule), 
    // and the actual foot pos (stance: the current pos, swing: the foot hold) at each dt_mpc
    // In engineering, just compute B matrix once
    Matrix3d R_WB, I_WB_inv;
    R_WB = AngleAxisd(euler_WB(2), Vector3d::UnitZ()) *
            AngleAxisd(euler_WB(1), Vector3d::UnitY()) *
            AngleAxisd(euler_WB(0), Vector3d::UnitX());
    I_WB_inv = (R_WB * _I_BB * R_WB.transpose()).inverse();
    for (int i = 0; i < 4; i++)
    {
        _Bc_mat.block<3,3>(6,3*i) = I_WB_inv * skew(R_WB * p_BF_list[i]);
        _Bc_mat.block<3,3>(9,3*i) = (1./_m) * Matrix3d::Identity();
        // TODO
        // stance: current foot pos
        // if (contact_schedule[i]) 
        // {
            
        // }
        // // swing: next foot hold 
        // else
        // {

        // }
    }
    // discretization
    _Bd_mat = _dt_mpc * _Bc_mat;
}

void ConvexMpcSolver::_calc_horizon_state(VectorXd body_state, VectorXd body_state_d)
{
    // update the current state (rpy pos w v g)
    _cur_state.segment<12>(0) = body_state;
    _cur_state(12) = -9.81;
    
    // update the desired mpc state trajectory
    for (int i = 0; i < HORIZON; i++)
    {
        _mpc_state_d.segment<2>(STATE_DIM*i) = body_state_d.segment<2>(0);
        _mpc_state_d(STATE_DIM*i + 2) = body_state_d(2) + _dt_mpc * i * body_state_d(8);
        _mpc_state_d.segment<2>(STATE_DIM*i + 3) = body_state_d.segment<2>(3) + 
                                                    _dt_mpc * i * body_state_d.segment<2>(9);
        _mpc_state_d(STATE_DIM*i + 5) = body_state_d(5);
        _mpc_state_d.segment<6>(STATE_DIM*i + 6) = body_state_d.segment<6>(6);
        _mpc_state_d(STATE_DIM*i + 12) = -9.81;
    }
}

void ConvexMpcSolver::_calc_qp_mat(vector<bool> contact_schedule)
{
    /*
    qp fromulation
    min  1/2 * U' * H * U + U' * g
    s.t. lb <= A * U <= ub    

    A_qp = [A0
            A1*A0
            A2*A1*A0
            ...
            An-1*An-2*...*A0]
    B_qp = [B0, 
            A1*B0,              B1, 
            A2*A1*B0,           A2*B1,              B2,
            ...
            An-1*An-2*...A1*B0, An-1*An-2*...A2*B1, ... , Bn-1]

    approximately 
    A_qp = [A
            A^2
            A^3
            ...
            A^n] 
    B_qp = [B, 
            A*B,     B, 
            A^2*B,   A*B,     B,
            ...
            A^n-1*B, A^n-2*B, ... , B]    
    */
    MatrixXd dense_hessian, dense_linear_constraint;

    for (int i = 0; i < HORIZON; i++){
        if (i == 0){
            _Aqp_mat.block<STATE_DIM, STATE_DIM>(0,0) = _Ad_mat;
        }
        else{
            _Aqp_mat.block<STATE_DIM, STATE_DIM>(STATE_DIM*i, 0) = _Ad_mat * 
                                _Aqp_mat.block<STATE_DIM, STATE_DIM>(STATE_DIM*(i-1), 0);
        }
        for (int j = 0; j < i+1; j++){
            if (j == i){
                _Bqp_mat.block<STATE_DIM, NUM_DOF>(STATE_DIM*i, NUM_DOF*i) = _Bd_mat;
            }
            else{
                _Bqp_mat.block<STATE_DIM, NUM_DOF>(STATE_DIM*i, NUM_DOF*j) = _Ad_mat *
                                _Bqp_mat.block<STATE_DIM, NUM_DOF>(STATE_DIM*(i-1), NUM_DOF*j);
            }
        }
    }
    dense_hessian = _Bqp_mat.transpose() * _Q_mat * _Bqp_mat;
    dense_hessian += _R_mat;
    _hessian_mat = dense_hessian.sparseView();
    _gradient_mat = _Bqp_mat.transpose() * _Q_mat * (_Aqp_mat * _cur_state - _mpc_state_d);

    // Theoretically, "horizon" contact schedule state should be computed
    // In eningeering, do not update the contact schedule
    dense_linear_constraint.resize(CONSTRAINT_DIM * HORIZON, NUM_DOF * HORIZON);
    for (int i = 0; i < HORIZON; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            // TODO: add slope
            dense_linear_constraint.block<5,3>(CONSTRAINT_DIM*i + 5*j, NUM_DOF*i + 3*j) = _linear_constraint_leg;
        }
    }
    _linear_constraint_mat = dense_linear_constraint.sparseView();
    
    for (int i = 0; i < HORIZON; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            if (contact_schedule[j]) // stance
            {
                _lb_mat(CONSTRAINT_DIM*i + 5*j + 4) = _Fz_min;
                _ub_mat(CONSTRAINT_DIM*i + 5*j + 4) = _Fz_max;
            }
            else // swing
            {
                _ub_mat(CONSTRAINT_DIM*i + 5*j + 4) = 0.;
            }
        }
    }
}

VectorXd ConvexMpcSolver::solve(VectorXd body_state, 
                                VectorXd body_state_d, 
                                vector<Vector3d> p_BF_list, 
                                vector<bool> contact_schedule)
{    
    reset();
    _calc_A_mat(body_state(2));
    _calc_B_mat(body_state.segment<3>(0), p_BF_list, contact_schedule);
    _calc_horizon_state(body_state, body_state_d);
    _calc_qp_mat(contact_schedule);

    // logger->debug("cur_state:\n{}", _cur_state.transpose());
    // logger->debug("mpc_state:\n{}", _mpc_state_d.transpose());
    // logger->debug("Aqp:\n{}", _Aqp_mat);
    // logger->debug("Bqp:\n{}", _Bqp_mat);
    // logger->debug("linear_constraint:\n{}", _linear_constraint_mat.toDense());
    // logger->debug("lb:\n{}", _lb_mat.transpose());
    // logger->debug("ub:\n{}", _ub_mat.transpose());

    // construct solver
    OsqpEigen::Solver qp_solver;
    // settings
    qp_solver.settings()->setVerbosity(false); // print the details of solution 
    qp_solver.settings()->setWarmStart(false); 
    // set data
    qp_solver.data()->setNumberOfVariables(NUM_DOF * HORIZON); // 120
    qp_solver.data()->setNumberOfConstraints(CONSTRAINT_DIM * HORIZON); // 200
    qp_solver.data()->setHessianMatrix(_hessian_mat);
    qp_solver.data()->setGradient(_gradient_mat);
    qp_solver.data()->setLinearConstraintsMatrix(_linear_constraint_mat);
    qp_solver.data()->setLowerBound(_lb_mat);
    qp_solver.data()->setUpperBound(_ub_mat);
    // init solver
    if (!qp_solver.initSolver())
    {
        return VectorXd::Zero(NUM_DOF,1);
    }
    // solve
    if (qp_solver.solveProblem() != OsqpEigen::ErrorExitFlag::NoError)
    {
        return VectorXd::Zero(NUM_DOF,1);
    }

    return qp_solver.getSolution().block<NUM_DOF,1>(0,0);
}
