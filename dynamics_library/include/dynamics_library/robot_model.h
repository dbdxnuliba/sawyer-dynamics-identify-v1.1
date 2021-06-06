#ifndef ROBOT_MODEL_H_
#define ROBOT_MODEL_H_

#include <Eigen/Dense>

using namespace Eigen;

namespace robot_dyn
{
class RobotModel
{
public:
    RobotModel(){};
    ~RobotModel(){};

    void InitModel();

    // Mi MXi MYi MZi XXi XYi XZi YYi YZi ZZi FVi FSi
    void SetDynamicsParameters(const VectorXd param);

    VectorXd calcu_inv_dyn(const VectorXd q, const VectorXd qDot, const VectorXd qDDot);
    MatrixXd calcu_InertiaMatrix(const VectorXd q);
    VectorXd calcu_CoriolisCentripetal(const VectorXd q, const VectorXd qDot);
    VectorXd calcu_Gravity(const VectorXd q);

    int dof;
    double g;

    VectorXd qMin; VectorXd qMax; 
    VectorXd qDotMin; VectorXd qDotMax; 
    VectorXd qDDotMin; VectorXd qDDotMax;

    // Mi MXi MYi MZi XXi XYi XZi YYi YZi ZZi FVi FSi
    unsigned int Psi_num;
    unsigned int Ps_num;
    VectorXi Ps_flag;
    VectorXd Ps;
    unsigned int Pb_num;
    VectorXd Pb;
    MatrixXd R1;
    MatrixXd R2;
    double qr_threshold;

private:
    double M1; double MX1; double MY1; double MZ1; double XX1; double XY1; double XZ1; double YY1; double YZ1; double ZZ1; double FV1; double FS1;
    double M2; double MX2; double MY2; double MZ2; double XX2; double XY2; double XZ2; double YY2; double YZ2; double ZZ2; double FV2; double FS2;
    double M3; double MX3; double MY3; double MZ3; double XX3; double XY3; double XZ3; double YY3; double YZ3; double ZZ3; double FV3; double FS3;
    double M4; double MX4; double MY4; double MZ4; double XX4; double XY4; double XZ4; double YY4; double YZ4; double ZZ4; double FV4; double FS4;
    double M5; double MX5; double MY5; double MZ5; double XX5; double XY5; double XZ5; double YY5; double YZ5; double ZZ5; double FV5; double FS5;
    double M6; double MX6; double MY6; double MZ6; double XX6; double XY6; double XZ6; double YY6; double YZ6; double ZZ6; double FV6; double FS6;
    double M7; double MX7; double MY7; double MZ7; double XX7; double XY7; double XZ7; double YY7; double YZ7; double ZZ7; double FV7; double FS7;

    VectorXd tau;
};

}

#endif
