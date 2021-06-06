#include <dynamics_library/robot_model.h>
#include <Eigen/Dense>
#include <iostream>
#include <math.h>
#include <time.h>

using namespace Eigen;

int main(int argc, char ** argv)
{
    robot_dyn::RobotModel robot;
    robot.InitModel();

    VectorXd param = VectorXd::Zero(robot.Ps_num);
    param << 1, 0, 0, 0, 0, 0, 0, 0, 0, 0.269807, 10.0627, 0.000811057, 1, 0.0142309, 1.00035, 0, 0.301216, -0.00931855, 0.0135568, 0, 0.0296092, 0.3316, 5.09696, -0.036601, 1, 0.00781326, 0.255642, 0, 0.0172233, 0.000531151, 0.000958207, 0, 0.0196323, 0.0238513, 5.01698, 0.0134498, 1, 0.00774561, 0.108827, 0, -0.0226316, -0.00578521, -0.0139302, 0, -0.00307844, -0.0130958, 0.954849, -0.024785, 1, 0.00309203, -0.0138228, 0, 0.00935767, 0.00607404, -0.0034507, 0, -0.00659333, 0.0149614, 1.99832, 0.014847, 1, -0.0148274, 0.152681, 0, 0.0720264, -0.0105593, 0.00568894, 0, -0.00166185, 0.0727083, 1.0371, -0.00795955, 1, 0.000430787, 0.00482754, 0, -0.00727855, -0.00782562, -0.00779483, 0, 0.00748447, 0.0496921, 1.00892, -0.000381248;
    robot.SetDynamicsParameters(param);

    MatrixXd M = MatrixXd::Zero(robot.dof,robot.dof);
    VectorXd H = VectorXd::Zero(robot.dof);
    VectorXd G = VectorXd::Zero(robot.dof);

    clock_t start, finish;
    double totaltime;

    for (unsigned int i=0; i<100; i++)
    {
        VectorXd q = VectorXd::Random(robot.dof);
        VectorXd qDot = VectorXd::Random(robot.dof);

        start = clock();

        M = robot.calcu_InertiaMatrix(q);
        H = robot.calcu_CoriolisCentripetal(q, qDot);
        G = robot.calcu_Gravity(q);

        finish = clock();
        totaltime = (double)(finish-start)/CLOCKS_PER_SEC;
        std::cout << "totaltime: " << totaltime*1000.0 << " ms" << std::endl;
    }

    for (unsigned int i=0; i<10; i++)
    {
        VectorXd q = VectorXd::Random(robot.dof);
        VectorXd qDot = VectorXd::Random(robot.dof);

        M = robot.calcu_InertiaMatrix(q);
        H = robot.calcu_CoriolisCentripetal(q, qDot);
        G = robot.calcu_Gravity(q);

        std::cout << "M:" << std::endl << M << std::endl;
        std::cout << "H:" << std::endl << H << std::endl;
        std::cout << "G:" << std::endl << G << std::endl;
    }

    return 0;
}
