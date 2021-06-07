#include <dynamics_library/robot_model.h>
#include <dynamics_library/robot_math.h>
#include <math.h>

namespace robot_dyn
{
void RobotModel::InitModel()
{
    dof = 7;
    g = -9.81;

    qMin.resize(dof);
    qMax.resize(dof);
    qDotMin.resize(dof);
    qDotMax.resize(dof);
    qDDotMin.resize(dof);
    qDDotMax.resize(dof);

    Psi_num = 12;
    Ps_num = Psi_num*dof;
    Ps_flag.resize(Ps_num);
    Ps.resize(Ps_num);
    Pb_num = 0;
    qr_threshold = 1.0e-10;

    tau.resize(dof);
}

// Mi MXi MYi MZi XXi XYi XZi YYi YZi ZZi FVi FSi
void RobotModel::SetDynamicsParameters(const VectorXd param)
{
    M1 = param(0);           MX1 = param(1);           MY1 = param(2);           MZ1 = param(3);           XX1 = param(4);           XY1 = param(5);           XZ1 = param(6);           YY1 = param(7);           YZ1 = param(8);           ZZ1 = param(9);           FV1 = param(10);           FS1 = param(11);
    M2 = param(Psi_num+0);   MX2 = param(Psi_num+1);   MY2 = param(Psi_num+2);   MZ2 = param(Psi_num+3);   XX2 = param(Psi_num+4);   XY2 = param(Psi_num+5);   XZ2 = param(Psi_num+6);   YY2 = param(Psi_num+7);   YZ2 = param(Psi_num+8);   ZZ2 = param(Psi_num+9);   FV2 = param(Psi_num+10);   FS2 = param(Psi_num+11);
    M3 = param(2*Psi_num+0); MX3 = param(2*Psi_num+1); MY3 = param(2*Psi_num+2); MZ3 = param(2*Psi_num+3); XX3 = param(2*Psi_num+4); XY3 = param(2*Psi_num+5); XZ3 = param(2*Psi_num+6); YY3 = param(2*Psi_num+7); YZ3 = param(2*Psi_num+8); ZZ3 = param(2*Psi_num+9); FV3 = param(2*Psi_num+10); FS3 = param(2*Psi_num+11);
    M4 = param(3*Psi_num+0); MX4 = param(3*Psi_num+1); MY4 = param(3*Psi_num+2); MZ4 = param(3*Psi_num+3); XX4 = param(3*Psi_num+4); XY4 = param(3*Psi_num+5); XZ4 = param(3*Psi_num+6); YY4 = param(3*Psi_num+7); YZ4 = param(3*Psi_num+8); ZZ4 = param(3*Psi_num+9); FV4 = param(3*Psi_num+10); FS4 = param(3*Psi_num+11);
    M5 = param(4*Psi_num+0); MX5 = param(4*Psi_num+1); MY5 = param(4*Psi_num+2); MZ5 = param(4*Psi_num+3); XX5 = param(4*Psi_num+4); XY5 = param(4*Psi_num+5); XZ5 = param(4*Psi_num+6); YY5 = param(4*Psi_num+7); YZ5 = param(4*Psi_num+8); ZZ5 = param(4*Psi_num+9); FV5 = param(4*Psi_num+10); FS5 = param(4*Psi_num+11);
    M6 = param(5*Psi_num+0); MX6 = param(5*Psi_num+1); MY6 = param(5*Psi_num+2); MZ6 = param(5*Psi_num+3); XX6 = param(5*Psi_num+4); XY6 = param(5*Psi_num+5); XZ6 = param(5*Psi_num+6); YY6 = param(5*Psi_num+7); YZ6 = param(5*Psi_num+8); ZZ6 = param(5*Psi_num+9); FV6 = param(5*Psi_num+10); FS6 = param(5*Psi_num+11);
    M7 = param(6*Psi_num+0); MX7 = param(6*Psi_num+1); MY7 = param(6*Psi_num+2); MZ7 = param(6*Psi_num+3); XX7 = param(6*Psi_num+4); XY7 = param(6*Psi_num+5); XZ7 = param(6*Psi_num+6); YY7 = param(6*Psi_num+7); YZ7 = param(6*Psi_num+8); ZZ7 = param(6*Psi_num+9); FV7 = param(6*Psi_num+10); FS7 = param(6*Psi_num+11);
}

VectorXd RobotModel::calcu_inv_dyn(const VectorXd q, const VectorXd qDot, const VectorXd qDDot)
{
    double th1 = q(0);
    double th2 = q(1)-M_PI_2;
    double th3 = q(2);
    double th4 = q(3);
    double th5 = q(4);
    double th6 = q(5);
    double th7 = q(6);
    double QP1 = qDot(0);
    double QP2 = qDot(1);
    double QP3 = qDot(2);
    double QP4 = qDot(3);
    double QP5 = qDot(4);
    double QP6 = qDot(5);
    double QP7 = qDot(6);
    double QDP1 = qDDot(0);
    double QDP2 = qDDot(1);
    double QDP3 = qDDot(2);
    double QDP4 = qDDot(3);
    double QDP5 = qDDot(4);
    double QDP6 = qDDot(5);
    double QDP7 = qDDot(6);
    double GZ = g;

    double C1 = cos(th1);
    double S1 = sin(th1);
    double C2 = cos(th2);
    double S2 = sin(th2);
    double C3 = cos(th3);
    double S3 = sin(th3);
    double C4 = cos(th4);
    double S4 = sin(th4);
    double C5 = cos(th5);
    double S5 = sin(th5);
    double C6 = cos(th6);
    double S6 = sin(th6);
    double C7 = cos(th7);
    double S7 = sin(th7);
    double DV61 = pow(QP1,2);
    double W12 = -QP1*S2;
    double W22 = -C2*QP1;
    double WP12 = -QDP1*S2 + QP2*W22;
    double WP22 = -C2*QDP1 - QP2*W12;
    double DV12 = pow(W12,2);
    double DV22 = W12*W22;
    double DV32 = QP2*W12;
    double DV42 = pow(W22,2);
    double DV52 = QP2*W22;
    double DV62 = pow(QP2,2);
    double U112 = -DV42 - DV62;
    double U212 = DV22 + QDP2;
    double U312 = DV32 - WP22;
    double U122 = DV22 - QDP2;
    double U222 = -DV12 - DV62;
    double U322 = DV52 + WP12;
    double U132 = DV32 + WP22;
    double U232 = DV52 - WP12;
    double U332 = -DV12 - DV42;
    double VSP12 = -0.081*DV61 - 0.1925*QDP1;
    double VSP22 = -0.1925*DV61 + 0.081*QDP1;
    double VP12 = C2*VSP12 + GZ*S2;
    double VP22 = C2*GZ - S2*VSP12;
    double W13 = C3*W12 - QP2*S3;
    double W23 = -C3*QP2 - S3*W12;
    double W33 = QP3 + W22;
    double WP13 = C3*WP12 - QDP2*S3 + QP3*W23;
    double WP23 = -C3*QDP2 - QP3*W13 - S3*WP12;
    double WP33 = QDP3 + WP22;
    double DV13 = pow(W13,2);
    double DV23 = W13*W23;
    double DV33 = W13*W33;
    double DV43 = pow(W23,2);
    double DV53 = W23*W33;
    double DV63 = pow(W33,2);
    double U113 = -DV43 - DV63;
    double U213 = DV23 + WP33;
    double U313 = DV33 - WP23;
    double U123 = DV23 - WP33;
    double U223 = -DV13 - DV63;
    double U323 = DV53 + WP13;
    double U133 = DV33 + WP23;
    double U233 = DV53 - WP13;
    double U333 = -DV13 - DV43;
    double VSP13 = 0.4*U122 + VP12;
    double VSP23 = 0.4*U222 + VP22;
    double VSP33 = 0.4*U322 + VSP22;
    double VP13 = C3*VSP13 - S3*VSP33;
    double VP23 = -C3*VSP33 - S3*VSP13;
    double W14 = C4*W13 + S4*W33;
    double W24 = C4*W33 - S4*W13;
    double W34 = QP4 - W23;
    double WP14 = C4*WP13 + QP4*W24 + S4*WP33;
    double WP24 = C4*WP33 - QP4*W14 - S4*WP13;
    double WP34 = QDP4 - WP23;
    double DV14 = pow(W14,2);
    double DV24 = W14*W24;
    double DV34 = W14*W34;
    double DV44 = pow(W24,2);
    double DV54 = W24*W34;
    double DV64 = pow(W34,2);
    double U114 = -DV44 - DV64;
    double U214 = DV24 + WP34;
    double U314 = DV34 - WP24;
    double U124 = DV24 - WP34;
    double U224 = -DV14 - DV64;
    double U324 = DV54 + WP14;
    double U134 = DV34 + WP24;
    double U234 = DV54 - WP14;
    double U334 = -DV14 - DV44;
    double VSP14 = 0.1685*U123 + VP13;
    double VSP24 = 0.1685*U223 + VP23;
    double VSP34 = 0.1685*U323 + VSP23;
    double VP14 = C4*VSP14 + S4*VSP34;
    double VP24 = C4*VSP34 - S4*VSP14;
    double W15 = C5*W14 - S5*W34;
    double W25 = -C5*W34 - S5*W14;
    double W35 = QP5 + W24;
    double WP15 = C5*WP14 + QP5*W25 - S5*WP34;
    double WP25 = -C5*WP34 - QP5*W15 - S5*WP14;
    double WP35 = QDP5 + WP24;
    double DV15 = pow(W15,2);
    double DV25 = W15*W25;
    double DV35 = W15*W35;
    double DV45 = pow(W25,2);
    double DV55 = W25*W35;
    double DV65 = pow(W35,2);
    double U115 = -DV45 - DV65;
    double U215 = DV25 + WP35;
    double U315 = DV35 - WP25;
    double U125 = DV25 - WP35;
    double U225 = -DV15 - DV65;
    double U325 = DV55 + WP15;
    double U135 = DV35 + WP25;
    double U235 = DV55 - WP15;
    double U335 = -DV15 - DV45;
    double VSP15 = 0.4*U124 + VP14;
    double VSP25 = 0.4*U224 + VP24;
    double VSP35 = 0.4*U324 - VSP24;
    double VP15 = C5*VSP15 - S5*VSP35;
    double VP25 = -C5*VSP35 - S5*VSP15;
    double W16 = C6*W15 + S6*W35;
    double W26 = C6*W35 - S6*W15;
    double W36 = QP6 - W25;
    double WP16 = C6*WP15 + QP6*W26 + S6*WP35;
    double WP26 = C6*WP35 - QP6*W16 - S6*WP15;
    double WP36 = QDP6 - WP25;
    double DV16 = pow(W16,2);
    double DV26 = W16*W26;
    double DV36 = W16*W36;
    double DV46 = pow(W26,2);
    double DV56 = W26*W36;
    double DV66 = pow(W36,2);
    double U116 = -DV46 - DV66;
    double U216 = DV26 + WP36;
    double U316 = DV36 - WP26;
    double U126 = DV26 - WP36;
    double U226 = -DV16 - DV66;
    double U326 = DV56 + WP16;
    double U136 = DV36 + WP26;
    double U236 = DV56 - WP16;
    double U336 = -DV16 - DV46;
    double VSP16 = -0.1363*U125 + VP15;
    double VSP26 = -0.1363*U225 + VP25;
    double VSP36 = -0.1363*U325 + VSP25;
    double VP16 = C6*VSP16 + S6*VSP36;
    double VP26 = C6*VSP36 - S6*VSP16;
    double W17 = C7*W16 - S7*W36;
    double W27 = -C7*W36 - S7*W16;
    double W37 = QP7 + W26;
    double WP17 = C7*WP16 + QP7*W27 - S7*WP36;
    double WP27 = -C7*WP36 - QP7*W17 - S7*WP16;
    double WP37 = QDP7 + WP26;
    double DV17 = pow(W17,2);
    double DV27 = W17*W27;
    double DV37 = W17*W37;
    double DV47 = pow(W27,2);
    double DV57 = W27*W37;
    double DV67 = pow(W37,2);
    double U117 = -DV47 - DV67;
    double U217 = DV27 + WP37;
    double U317 = DV37 - WP27;
    double U127 = DV27 - WP37;
    double U227 = -DV17 - DV67;
    double U327 = DV57 + WP17;
    double U137 = DV37 + WP27;
    double U237 = DV57 - WP17;
    double U337 = -DV17 - DV47;
    double VP17 = C7*VP16 + S7*VSP26;
    double VP27 = C7*VSP26 - S7*VP16;
    double F11 = -DV61*MX1 - MY1*QDP1;
    double F21 = -DV61*MY1 + MX1*QDP1;
    double F31 = -GZ*M1;
    double PSI11 = QP1*XZ1;
    double PSI21 = QP1*YZ1;
    double PSI31 = QP1*ZZ1;
    double No11 = -PSI21*QP1 + QDP1*XZ1;
    double No21 = PSI11*QP1 + QDP1*YZ1;
    double No31 = QDP1*ZZ1;
    double F12 = M2*VP12 + MX2*U112 + MY2*U122 + MZ2*U132;
    double F22 = M2*VP22 + MX2*U212 + MY2*U222 + MZ2*U232;
    double F32 = M2*VSP22 + MX2*U312 + MY2*U322 + MZ2*U332;
    double PSI12 = QP2*XZ2 + W12*XX2 + W22*XY2;
    double PSI22 = QP2*YZ2 + W12*XY2 + W22*YY2;
    double PSI32 = QP2*ZZ2 + W12*XZ2 + W22*YZ2;
    double No12 = -PSI22*QP2 + PSI32*W22 + QDP2*XZ2 + WP12*XX2 + WP22*XY2;
    double No22 = PSI12*QP2 - PSI32*W12 + QDP2*YZ2 + WP12*XY2 + WP22*YY2;
    double No32 = -PSI12*W22 + PSI22*W12 + QDP2*ZZ2 + WP12*XZ2 + WP22*YZ2;
    double F13 = M3*VP13 + MX3*U113 + MY3*U123 + MZ3*U133;
    double F23 = M3*VP23 + MX3*U213 + MY3*U223 + MZ3*U233;
    double F33 = M3*VSP23 + MX3*U313 + MY3*U323 + MZ3*U333;
    double PSI13 = W13*XX3 + W23*XY3 + W33*XZ3;
    double PSI23 = W13*XY3 + W23*YY3 + W33*YZ3;
    double PSI33 = W13*XZ3 + W23*YZ3 + W33*ZZ3;
    double No13 = -PSI23*W33 + PSI33*W23 + WP13*XX3 + WP23*XY3 + WP33*XZ3;
    double No23 = PSI13*W33 - PSI33*W13 + WP13*XY3 + WP23*YY3 + WP33*YZ3;
    double No33 = -PSI13*W23 + PSI23*W13 + WP13*XZ3 + WP23*YZ3 + WP33*ZZ3;
    double F14 = M4*VP14 + MX4*U114 + MY4*U124 + MZ4*U134;
    double F24 = M4*VP24 + MX4*U214 + MY4*U224 + MZ4*U234;
    double F34 = -M4*VSP24 + MX4*U314 + MY4*U324 + MZ4*U334;
    double PSI14 = W14*XX4 + W24*XY4 + W34*XZ4;
    double PSI24 = W14*XY4 + W24*YY4 + W34*YZ4;
    double PSI34 = W14*XZ4 + W24*YZ4 + W34*ZZ4;
    double No14 = -PSI24*W34 + PSI34*W24 + WP14*XX4 + WP24*XY4 + WP34*XZ4;
    double No24 = PSI14*W34 - PSI34*W14 + WP14*XY4 + WP24*YY4 + WP34*YZ4;
    double No34 = -PSI14*W24 + PSI24*W14 + WP14*XZ4 + WP24*YZ4 + WP34*ZZ4;
    double F15 = M5*VP15 + MX5*U115 + MY5*U125 + MZ5*U135;
    double F25 = M5*VP25 + MX5*U215 + MY5*U225 + MZ5*U235;
    double F35 = M5*VSP25 + MX5*U315 + MY5*U325 + MZ5*U335;
    double PSI15 = W15*XX5 + W25*XY5 + W35*XZ5;
    double PSI25 = W15*XY5 + W25*YY5 + W35*YZ5;
    double PSI35 = W15*XZ5 + W25*YZ5 + W35*ZZ5;
    double No15 = -PSI25*W35 + PSI35*W25 + WP15*XX5 + WP25*XY5 + WP35*XZ5;
    double No25 = PSI15*W35 - PSI35*W15 + WP15*XY5 + WP25*YY5 + WP35*YZ5;
    double No35 = -PSI15*W25 + PSI25*W15 + WP15*XZ5 + WP25*YZ5 + WP35*ZZ5;
    double F16 = M6*VP16 + MX6*U116 + MY6*U126 + MZ6*U136;
    double F26 = M6*VP26 + MX6*U216 + MY6*U226 + MZ6*U236;
    double F36 = -M6*VSP26 + MX6*U316 + MY6*U326 + MZ6*U336;
    double PSI16 = W16*XX6 + W26*XY6 + W36*XZ6;
    double PSI26 = W16*XY6 + W26*YY6 + W36*YZ6;
    double PSI36 = W16*XZ6 + W26*YZ6 + W36*ZZ6;
    double No16 = -PSI26*W36 + PSI36*W26 + WP16*XX6 + WP26*XY6 + WP36*XZ6;
    double No26 = PSI16*W36 - PSI36*W16 + WP16*XY6 + WP26*YY6 + WP36*YZ6;
    double No36 = -PSI16*W26 + PSI26*W16 + WP16*XZ6 + WP26*YZ6 + WP36*ZZ6;
    double F17 = M7*VP17 + MX7*U117 + MY7*U127 + MZ7*U137;
    double F27 = M7*VP27 + MX7*U217 + MY7*U227 + MZ7*U237;
    double F37 = M7*VP26 + MX7*U317 + MY7*U327 + MZ7*U337;
    double PSI17 = W17*XX7 + W27*XY7 + W37*XZ7;
    double PSI27 = W17*XY7 + W27*YY7 + W37*YZ7;
    double PSI37 = W17*XZ7 + W27*YZ7 + W37*ZZ7;
    double No17 = -PSI27*W37 + PSI37*W27 + WP17*XX7 + WP27*XY7 + WP37*XZ7;
    double No27 = PSI17*W37 - PSI37*W17 + WP17*XY7 + WP27*YY7 + WP37*YZ7;
    double No37 = -PSI17*W27 + PSI27*W17 + WP17*XZ7 + WP27*YZ7 + WP37*ZZ7;
    double N17 = MY7*VP26 - MZ7*VP27 + No17;
    double N27 = -MX7*VP26 + MZ7*VP17 + No27;
    double N37 = MX7*VP27 - MY7*VP17 + No37;
    double FDI17 = C7*F17 - F27*S7;
    double FDI37 = -C7*F27 - F17*S7;
    double E16 = F16 + FDI17;
    double E26 = F26 + F37;
    double E36 = F36 + FDI37;
    double N16 = C7*N17 - MY6*VSP26 - MZ6*VP26 - N27*S7 + No16;
    double N26 = MX6*VSP26 + MZ6*VP16 + N37 + No26;
    double N36 = -C7*N27 + MX6*VP26 - MY6*VP16 - N17*S7 + No36;
    double FDI16 = C6*E16 - E26*S6;
    double FDI36 = C6*E26 + E16*S6;
    double E15 = F15 + FDI16;
    double E25 = -E36 + F25;
    double E35 = F35 + FDI36;
    double N15 = C6*N16 - 0.1363*FDI36 + MY5*VSP25 - MZ5*VP25 - N26*S6 + No15;
    double N25 = -MX5*VSP25 + MZ5*VP15 - N36 + No25;
    double N35 = C6*N26 + 0.1363*FDI16 + MX5*VP25 - MY5*VP15 + N16*S6 + No35;
    double FDI15 = C5*E15 - E25*S5;
    double FDI35 = -C5*E25 - E15*S5;
    double E14 = F14 + FDI15;
    double E24 = E35 + F24;
    double E34 = F34 + FDI35;
    double N14 = C5*N15 + 0.4*FDI35 - MY4*VSP24 - MZ4*VP24 - N25*S5 + No14;
    double N24 = MX4*VSP24 + MZ4*VP14 + N35 + No24;
    double N34 = -C5*N25 - 0.4*FDI15 + MX4*VP24 - MY4*VP14 - N15*S5 + No34;
    double FDI14 = C4*E14 - E24*S4;
    double FDI34 = C4*E24 + E14*S4;
    double E13 = F13 + FDI14;
    double E23 = -E34 + F23;
    double E33 = F33 + FDI34;
    double N13 = C4*N14 + 0.1685*FDI34 + MY3*VSP23 - MZ3*VP23 - N24*S4 + No13;
    double N23 = -MX3*VSP23 + MZ3*VP13 - N34 + No23;
    double N33 = C4*N24 - 0.1685*FDI14 + MX3*VP23 - MY3*VP13 + N14*S4 + No33;
    double FDI13 = C3*E13 - E23*S3;
    double FDI33 = -C3*E23 - E13*S3;
    double E12 = F12 + FDI13;
    double E22 = E33 + F22;
    double E32 = F32 + FDI33;
    double N12 = C3*N13 + 0.4*FDI33 + MY2*VSP22 - MZ2*VP22 - N23*S3 + No12;
    double N22 = -MX2*VSP22 + MZ2*VP12 + N33 + No22;
    double N32 = -C3*N23 - 0.4*FDI13 + MX2*VP22 - MY2*VP12 - N13*S3 + No32;
    double FDI12 = C2*E12 - E22*S2;
    double FDI32 = -C2*E22 - E12*S2;
    double E11 = F11 + FDI12;
    double E21 = E32 + F21;
    double E31 = F31 + FDI32;
    double N11 = C2*N12 + 0.1925*FDI32 - GZ*MY1 - N22*S2 + No11;
    double N21 = -0.081*FDI32 + GZ*MX1 + N32 + No21;
    double N31 = -C2*N22 + 0.081*E32 - 0.1925*FDI12 - N12*S2 + No31;
    double FDI11 = C1*E11 - E21*S1;
    double FDI21 = C1*E21 + E11*S1;
    double GAM1 = FS1*sign(QP1) + FV1*QP1 + N31;
    double GAM2 = FS2*sign(QP2) + FV2*QP2 + N32;
    double GAM3 = FS3*sign(QP3) + FV3*QP3 + N33;
    double GAM4 = FS4*sign(QP4) + FV4*QP4 + N34;
    double GAM5 = FS5*sign(QP5) + FV5*QP5 + N35;
    double GAM6 = FS6*sign(QP6) + FV6*QP6 + N36;
    double GAM7 = FS7*sign(QP7) + FV7*QP7 + N37;

    tau << GAM1, GAM2, GAM3, GAM4, GAM5, GAM6, GAM7;

    return tau;
}

MatrixXd RobotModel::calcu_InertiaMatrix(const VectorXd q)
{
    VectorXd qDot = VectorXd::Zero(dof);

    double g_ = g;
    g = 0;

    double FS1_ = FS1;
    double FS2_ = FS2;
    double FS3_ = FS3;
    double FS4_ = FS4;
    double FS5_ = FS5;
    double FS6_ = FS6;
    double FS7_ = FS7;
    FS1 = 0;
    FS2 = 0;
    FS3 = 0;
    FS4 = 0;
    FS5 = 0;
    FS6 = 0;
    FS7 = 0;

    MatrixXd M = MatrixXd::Zero(dof,dof);

    for (unsigned int i=0; i< dof; i++)
    {
        VectorXd qDDot = VectorXd::Zero(dof);

        qDDot(i) = 1;
        M.col(i) = calcu_inv_dyn(q,qDot,qDDot);
    }

    g = g_;
    
    FS1 = FS1_;
    FS2 = FS2_;
    FS3 = FS3_;
    FS4 = FS4_;
    FS5 = FS5_;
    FS6 = FS6_;
    FS7 = FS7_;
    
    return M;
}

VectorXd RobotModel::calcu_CoriolisCentripetal(const VectorXd q, const VectorXd qDot)
{
    VectorXd qDDot = VectorXd::Zero(dof);

    double g_ = g;
    g = 0;

    double FV1_ = FV1;
    double FV2_ = FV2;
    double FV3_ = FV3;
    double FV4_ = FV4;
    double FV5_ = FV5;
    double FV6_ = FV6;
    double FV7_ = FV7;
    FV1 = 0;
    FV2 = 0;
    FV3 = 0;
    FV4 = 0;
    FV5 = 0;
    FV6 = 0;
    FV7 = 0;
    
    double FS1_ = FS1;
    double FS2_ = FS2;
    double FS3_ = FS3;
    double FS4_ = FS4;
    double FS5_ = FS5;
    double FS6_ = FS6;
    double FS7_ = FS7;
    FS1 = 0;
    FS2 = 0;
    FS3 = 0;
    FS4 = 0;
    FS5 = 0;
    FS6 = 0;
    FS7 = 0;

    VectorXd H = VectorXd::Zero(dof);

    H = calcu_inv_dyn(q,qDot,qDDot);

    g = g_;

    FV1 = FV1_;
    FV2 = FV2_;
    FV3 = FV3_;
    FV4 = FV4_;
    FV5 = FV5_;
    FV6 = FV6_;
    FV7 = FV7_;
    
    FS1 = FS1_;
    FS2 = FS2_;
    FS3 = FS3_;
    FS4 = FS4_;
    FS5 = FS5_;
    FS6 = FS6_;
    FS7 = FS7_;

    return H;
}

VectorXd RobotModel::calcu_Gravity(const VectorXd q)
{
    VectorXd qDot = VectorXd::Zero(dof);
    VectorXd qDDot = VectorXd::Zero(dof);

    double FS1_ = FS1;
    double FS2_ = FS2;
    double FS3_ = FS3;
    double FS4_ = FS4;
    double FS5_ = FS5;
    double FS6_ = FS6;
    double FS7_ = FS7;
    FS1 = 0;
    FS2 = 0;
    FS3 = 0;
    FS4 = 0;
    FS5 = 0;
    FS6 = 0;
    FS7 = 0;

    VectorXd G = VectorXd::Zero(dof);

    G = calcu_inv_dyn(q,qDot,qDDot);

    FS1 = FS1_;
    FS2 = FS2_;
    FS3 = FS3_;
    FS4 = FS4_;
    FS5 = FS5_;
    FS6 = FS6_;
    FS7 = FS7_;

    return G;
}

VectorXd RobotModel::calcu_Friction(const VectorXd qDot)
{
    VectorXd F = VectorXd::Zero(dof);

    F(0) = FV1*qDot(0)+FS1*sign(qDot(0));
    F(1) = FV2*qDot(1)+FS2*sign(qDot(1));
    F(2) = FV3*qDot(2)+FS3*sign(qDot(2));
    F(3) = FV4*qDot(3)+FS4*sign(qDot(3));
    F(4) = FV5*qDot(4)+FS5*sign(qDot(4));
    F(5) = FV6*qDot(5)+FS6*sign(qDot(5));
    F(6) = FV7*qDot(6)+FS7*sign(qDot(6));

    return F;
}

}









