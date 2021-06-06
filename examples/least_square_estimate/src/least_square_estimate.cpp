#include <least_square_estimate/deal_txt_data.h>
#include <string>
#include <Eigen/Dense>
#include <iostream>
#include <dynamics_library/robot_model.h>
#include <dynamics_library/robot_identify.h>
#include <fstream>

using namespace Eigen;

int main(int argc, char **argv)  
{  
    std::string file_path;
    unsigned int row;
    unsigned int col;

    file_path = "/home/stz/sawyer-dynamics-identify-v1.1/examples/least_square_estimate/filter_data/q_filter.txt";
    deal_txt::get_matrix_size(file_path, row, col);
    MatrixXd q_filter = MatrixXd::Zero(row,col);
    deal_txt::get_data_matrix(file_path, col, q_filter);

    file_path = "/home/stz/sawyer-dynamics-identify-v1.1/examples/least_square_estimate/filter_data/qDot_filter.txt";
    deal_txt::get_matrix_size(file_path, row, col);
    MatrixXd qDot_filter = MatrixXd::Zero(row,col);
    deal_txt::get_data_matrix(file_path, col, qDot_filter);

    file_path = "/home/stz/sawyer-dynamics-identify-v1.1/examples/least_square_estimate/filter_data/qDDot_filter.txt";
    deal_txt::get_matrix_size(file_path, row, col);
    MatrixXd qDDot_filter = MatrixXd::Zero(row,col);
    deal_txt::get_data_matrix(file_path, col, qDDot_filter);

    file_path = "/home/stz/sawyer-dynamics-identify-v1.1/examples/least_square_estimate/filter_data/tau_filter.txt";
    deal_txt::get_matrix_size(file_path, row, col);
    MatrixXd tau_filter = MatrixXd::Zero(row,col);
    deal_txt::get_data_matrix(file_path, col, tau_filter);

    robot_dyn::RobotModel robot;
    robot.InitModel();

    MatrixXd W = MatrixXd::Zero(row*robot.dof, robot.Ps_num);
    VectorXd tau = VectorXd::Zero(row*robot.dof);
    for (unsigned int i=0; i<row; i++)
    {
        W.middleRows(i*robot.dof, robot.dof) =
            robot_dyn::calcu_Ys(&robot, q_filter.row(i).transpose(), 
                qDot_filter.row(i).transpose(), qDDot_filter.row(i).transpose());

        tau.segment(i*robot.dof, robot.dof) = tau_filter.row(i).transpose();
    }

    MatrixXd Wb;
    qr_decompose(&robot, W, Wb);

    robot.Pb.resize(robot.Pb_num);
    robot.Pb = (Wb.transpose()*Wb).inverse()*Wb.transpose()*tau;
    
    // Pb = P1+inv(R1)*R2*P2
    VectorXd P1 = VectorXd::Zero(robot.Pb_num);
    VectorXd P2 = VectorXd::Zero(robot.Ps_num-robot.Pb_num);

    double m = 1.0;
    int P2_count = -1;
    for (unsigned int i=0; i< robot.Ps_num; i++)
    {
        if (robot.Ps_flag(i)==0)
        {
            P2_count += 1;
            if (i%robot.Psi_num==0)
            {
                P2(P2_count) = m;
            }
        }
    }

    P1 = robot.Pb-robot.R1.inverse()*robot.R2*P2;

    int P1_count = -1;
    P2_count = -1;
    for (unsigned int i=0; i< robot.Ps_num; i++)
    {
        if (robot.Ps_flag(i)==1)
        {
            P1_count += 1;
            robot.Ps(i) = P1(P1_count);
        }
        else
        {
            P2_count += 1;
            robot.Ps(i) = P2(P2_count);
        }
    }

    robot.SetDynamicsParameters(robot.Ps);

    MatrixXd tau_iden = MatrixXd::Zero(row, col);
    for (unsigned int i=0; i<row; i++)
    {
        tau_iden.row(i) = robot.calcu_inv_dyn(q_filter.row(i).transpose(), 
            qDot_filter.row(i).transpose(), qDDot_filter.row(i).transpose()).transpose();
    }

    std::ofstream dyn_params_file;
    dyn_params_file.open("/home/stz/sawyer-dynamics-identify-v1.1/examples/least_square_estimate/identify_result/dyn_params.txt");
    dyn_params_file << "Pb = [";
    for (unsigned int j=0; j<robot.Pb_num; j++)
    {
        dyn_params_file << robot.Pb(j);
        if (j!=(robot.Pb_num-1))
        {
            dyn_params_file << ", ";
        }
    }
    dyn_params_file << "]" << std::endl;
    dyn_params_file << "Ps = [";
    for (unsigned int j=0; j<robot.Ps_num; j++)
    {
        dyn_params_file << robot.Ps(j);
        if (j!=(robot.Ps_num-1))
        {
            dyn_params_file << ", ";
        }
    }
    dyn_params_file << "]" << std::endl;
    dyn_params_file.close();

    std::ofstream tau_iden_file;
    tau_iden_file.open("/home/stz/sawyer-dynamics-identify-v1.1/examples/least_square_estimate/identify_result/tau_iden.txt");
    tau_iden_file << "[";
    for (unsigned int i=0; i<row; i++)
    {
        for (unsigned int j=0; j<col; j++)
        {
            tau_iden_file << tau_iden(i,j);
            if (j!=(col-1))
            {
                tau_iden_file << ", ";
            }
            else
            {
                if (i!=(row-1))
                {
                    tau_iden_file << ";...";
                }
                else
                {
                    tau_iden_file << "]";
                }
            }
        }
        
        if (i!=(row-1))
        {
            tau_iden_file << std::endl;
        }
    }
    
    tau_iden_file.close();

    return 0;  
}  