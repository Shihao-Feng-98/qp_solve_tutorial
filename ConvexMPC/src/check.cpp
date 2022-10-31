#include <pthread.h> // -lpthread
#include <sys/mman.h> // mlockall(MCL_CURRENT|MCL_FUTURE)
#include <unistd.h> // sleep
using namespace std;
#include <periodic_rt_task.h>
#include <C_timer.h>
#include "ConvexMpcSolver.h"

double dt_mpc, R;
Matrix<double,13,1> Q;
Matrix<double, 12, 1> body_state, body_state_d; 
Vector3d p_BF_FR, p_BF_FL, p_BF_RR, p_BF_RL;
vector<Vector3d> p_BF_list;
vector<bool> contact_schedule;

void* main_loop(void* argc)
{
    CTimer timer_total;
    ConvexMpcSolver solver(dt_mpc, Q, R);
    int N = 10;

    cout << "[Main Thread]: thread start\n";

    timer_total.reset();
    for (int i = 0; i < N; i++)
    {
        solver.solve(body_state, body_state_d, p_BF_list, contact_schedule);
    }
    double duration = timer_total.end()/1000;

    cout << "Total time: " << duration << " ms\n";
    cout << "Solve time: " << duration/N << " ms\n"; //  一次7ms? 多次75ms?
    
    cout << "[Main Thread]: thread end\n";

    return nullptr;
}

int main(int argc, char **argv)
{
    dt_mpc = 0.002;
    Q << 25.,25.,10.,
        1., 1., 100.,
        0., 0., 0.3,
        0.2, 0.2, 20.,
        0.;
    R = 1e-5;
    body_state << 0.1, 0.1, 0.,
                0., 0., 0.23,
                0., 0., 1.,
                1., 1., 0.;
    body_state_d << 0., 0., 0.,
                    0., 0., 0.23,
                    0., 0., 1.,
                    1., 1., 0.;
    p_BF_FR << 0.2, -0.1, -0.23;
    p_BF_FL << 0.2, 0.1, -0.23;
    p_BF_RR << -0.2, -0.1, -0.23;
    p_BF_RL << -0.2, 0.1, -0.23;
    p_BF_list.emplace_back(p_BF_FR);
    p_BF_list.emplace_back(p_BF_FL);
    p_BF_list.emplace_back(p_BF_RR);
    p_BF_list.emplace_back(p_BF_RL);
    contact_schedule = {1,1,1,1};

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
