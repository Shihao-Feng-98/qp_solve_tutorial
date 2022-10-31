#include <balance_qp_solver.h>

BalanceQpControllerSolver::BalanceQpControllerSolver(double m, Vector3d I, 
                                                    double mu, double Fz_min, double Fz_max, 
                                                    VectorXd Q, double R)
{
    this->_m = m;
    this->_I_BB.diagonal() << I(0), I(1), I(2);
    this->_mu = mu;
    this->_Fz_min = Fz_min; 
    this->_Fz_max = Fz_max;
    this->_Q_mat.diagonal() << Q(0), Q(1), Q(2), Q(3), Q(4), Q(5);
    this->_R_mat = R * MatrixXd::Identity(12,12);

    this->_hessian_mat.resize(12,12);
    this->_gradient_mat.resize(12);
    this->_linear_mat.resize(20,12);
    // 这部分只能insert一次，否则会报错
    for (int i = 0; i < 4; i++) 
    {
        // 0 <= F_xi + u F_zi <= +oo
        _linear_mat.insert(5*i, 3*i) = 1;
        _linear_mat.insert(5*i, 3*i+2) = _mu;
        // 0 <= -F_xi + u F_zi <= +oo
        _linear_mat.insert(5*i+1, 3*i) = -1;
        _linear_mat.insert(5*i+1, 3*i+2) = _mu;
        // 0 <= F_yi + u F_zi <= +oo
        _linear_mat.insert(5*i+2, 3*i+1) = 1;
        _linear_mat.insert(5*i+2, 3*i+2) = _mu;
        // 0 <= -F_yi + u F_zi <= +oo
        _linear_mat.insert(5*i+3, 3*i+1) = -1;
        _linear_mat.insert(5*i+3, 3*i+2) = _mu;
        // F_zmin <= F_zi <= F_zmax
        _linear_mat.insert(5*i+4, 3*i+2) = 1;
    }
    this->_lb_mat.resize(20);
    this->_lb_mat.setZero();
    this->_ub_mat = OsqpEigen::INFTY * VectorXd::Ones(20);
}

BalanceQpControllerSolver::~BalanceQpControllerSolver()
{
    // nothing to do
}

VectorXd BalanceQpControllerSolver::solve(bool const *scheduled_contact, 
                                        const MatrixXd &A_mat, 
                                        const VectorXd &b_mat)
{
    // construct matrix
    MatrixXd dense_hessian = A_mat.transpose() * _Q_mat * A_mat + _R_mat;
    _hessian_mat = dense_hessian.sparseView();
    _gradient_mat = - A_mat.transpose() * _Q_mat * b_mat;
    // F_zmin <= F_zi <= F_zmax
    for (int i = 0; i < 4; i++)
    {
        double contact_flag = scheduled_contact[i] ? 1 : 0;
        _lb_mat(5*i+4) = contact_flag * _Fz_min;
        _ub_mat(5*i+4) = contact_flag * _Fz_max;
    }

    // construct solver
    OsqpEigen::Solver qp_solver;
    // settings
    qp_solver.settings()->setVerbosity(false); // 是否打印求解的详细结果
    qp_solver.settings()->setWarmStart(false); // ?
    // set data
    qp_solver.data()->setNumberOfVariables(12); // NUM_LEG * 3
    qp_solver.data()->setNumberOfConstraints(20); // NUM_LEG * 5
    qp_solver.data()->setHessianMatrix(_hessian_mat);
    qp_solver.data()->setGradient(_gradient_mat);
    qp_solver.data()->setLinearConstraintsMatrix(_linear_mat);
    qp_solver.data()->setLowerBound(_lb_mat);
    qp_solver.data()->setUpperBound(_ub_mat);
    // init solver
    if (!qp_solver.initSolver())
    {
        return VectorXd::Zero(12,1);
    }
    // solve
    if (qp_solver.solveProblem() != OsqpEigen::ErrorExitFlag::NoError)
    {
        return VectorXd::Zero(12,1);
    }

    return qp_solver.getSolution();
}
