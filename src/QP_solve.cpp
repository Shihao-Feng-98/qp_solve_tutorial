#include <iostream>
#include <vector>
#include <pthread.h> // -lpthread
#include <sys/mman.h> // mlockall(MCL_CURRENT|MCL_FUTURE)
#include <unistd.h> // sleep
using namespace std;
#include <Eigen/Dense>
using namespace Eigen;
#include "pinocchio/algorithm/kinematics.hpp"
#include "pinocchio/math/rpy.hpp"
namespace pin = pinocchio;
#include <balance_qp_solver.h>
#include <periodic_rt_task.h>
#include <C_timer.h>


double m, mu, Fz_min, Fz_max, R;
Vector3d I, g;
Matrix<double,6,1> Q;
DiagonalMatrix<double,3> I_BB;
DiagonalMatrix<double,3> kp_lin, kd_lin, kp_ang, kd_ang;
Vector3d ddp_WB_d, dp_WB_d, dp_WB, p_WB_d, p_WB;
Vector3d dw_WB_d, w_WB_d, w_WB, theta_WB_d, theta_WB;
vector<Vector3d> p_BH_vec(4,Vector3d::Zero()), p_HF_vec(4,Vector3d::Zero()), p_WF_vec(4,Vector3d::Zero());
Matrix3d R_WB;
bool scheduled_contact[4];
Matrix<double,6,12> A_mat;
Matrix<double,6,1> b_mat;


void calculation(BalanceQpControllerSolver &solver)
{
    // 实时更新的数据
    // 构造b矩阵
    // 身体PD控制
    p_WB_d << 0.001, 0., 0.25; // 0.3 * 0.0033 = 0.001
    dp_WB_d << 0.3, 0., 0.25;
    p_WB << 0., 0., 0.25;
    dp_WB << 0., 0., 0.25;

    theta_WB_d << 0., 0., 0.;
    w_WB_d << 0., 0., 0.;
    theta_WB << 0.05, 0., 0.; 
    w_WB << 0.3, 0., 0.;

    ddp_WB_d = kp_lin * (p_WB_d - p_WB) + kd_lin * (dp_WB_d - dp_WB);
    dw_WB_d = kp_ang * (theta_WB_d - theta_WB) + kd_ang * (w_WB_d - w_WB);
    
    b_mat.block<3,1>(0,0) = m * (ddp_WB_d - g);
    b_mat.block<3,1>(3,0) = I_BB * dw_WB_d;

    // 构造A矩阵
    // 足端力平衡
    // FR FL RR RL
    p_BH_vec[0] << 0.2155, -0.0625, 0.;
    p_BH_vec[1] << 0.2155, 0.0625, 0.;
    p_BH_vec[2] << -0.2155, -0.0625, 0.;
    p_BH_vec[3] << -0.2155, 0.0625, 0.;

    p_HF_vec[0] << 0., -0.15, -0.25;
    p_HF_vec[1] << 0., 0.15, -0.25;
    p_HF_vec[2] << 0., -0.15, -0.25;
    p_HF_vec[3] << 0., 0.15, -0.25;

    R_WB = pin::rpy::rpyToMatrix(theta_WB);
    for (int i = 0; i < 4; i++)
    {
        // assume R_WB * p_BF = p_WF
        p_WF_vec[i] = R_WB * (p_BH_vec[i] + p_HF_vec[i]);
    }

    // stand gait
    scheduled_contact[0] = 1;
    scheduled_contact[1] = 1;
    scheduled_contact[2] = 1;
    scheduled_contact[3] = 1;

    for (int i = 0; i < 4; i++)
    {
        A_mat.block<3,3>(0,3*i) = Matrix3d::Identity();
        A_mat.block<3,3>(3,3*i) = pin::skew(p_WF_vec[i]);
    }

    // 求解
    solver.solve(scheduled_contact, A_mat, b_mat);
}

void* main_loop(void* argc)
{
    CTimer timer_total;
    BalanceQpControllerSolver solver(m, I, mu, Fz_min, Fz_max, Q, R);
    int N = 10000;

    cout << "[Main Thread]: thread start\n";

    timer_total.reset();
    for (int i = 0; i < N; i++)
    {
        calculation(solver);
    }
    double duration = timer_total.end()/1000;

    cout << "Total time: " << duration << " ms\n";
    cout << "Solve time: " << duration/N << " ms\n"; // 0.05 ms以内
    
    cout << "[Main Thread]: thread end\n";

    return nullptr;
}

int main(int argc, char **argv)
{
    m = 12.6;
    mu = 0.6; 
    Fz_min = 1.;
    Fz_max = 150; 
    R = 1e-3;
    I << 0.133872, 0.450659, 0.469184;
    Q << 1., 1., 1., 100., 100., 25.;
    g << 0., 0., -9.81;
    I_BB.diagonal() << I(0), I(1), I(2);
    
    kp_lin.diagonal() << 10, 10, 100;
    kd_lin.diagonal() << 40, 40, 15;
    kp_ang.diagonal() << 48, 100, 1;
    kd_ang.diagonal() << 12, 18, 10;

    /*
    mlockall 锁定进程中所有映射到地址空间的页
    MCL_CURRENT 已经映射的进程地址，MCL_FUTURE 将来映射的进程地址
    */
    if(mlockall(MCL_CURRENT|MCL_FUTURE) == -1) {
        cout << "mlockall failed: %m\n"; // 输入上一个函数的错误信息
        return -2;
    }

    // 主控制线程
    PeriodicRtTask *main_task = new PeriodicRtTask("[Main Control Thread]", 95, main_loop);
    sleep(1); 
    // 析构函数会join线程，等待子线程结束
    delete main_task;

    return 0;
}
