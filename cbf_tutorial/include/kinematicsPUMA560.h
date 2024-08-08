#ifndef _KINEMATICSPUMA560_H_
#define _KINEMATICSPUMA560_H_

#include <functional>
#include <Eigen/Core>

#ifndef M_PI
   #define M_PI 3.14159265358979323846
#endif

// Constants for PUMA 560, assuming the following DH parameters (see "Introduction to robotics", Craig, 2005, page 80)
// The offset on q2 is consistent with the convention used in the original UNIMATION controllers (see footnote at page 78 of the same book)
/*
+-------+--------+---------+-------+-------+
| joint | d      | a       | alpha | offset|
+-------+--------+---------+-------+-------+
|    q0 | 0      | 0       | 0     | 0     |
|    q1 | 0      | 0       | -PI/2 | 0     |
|    q2 | 0.15   | 0.4318  | 0     | -PI   |
|    q3 | 0.4318 | 0.02    | -PI/2 | 0     |
|    q4 | 0      | 0       |  PI/2 | 0     |
|    q5 | 0      | 0       | -PI/2 | 0     |
+-------+--------+---------+-------+-------+
*/
const double d3 = 0.15;
const double d4 = 0.4318;
const double a2 = 0.4318;
const double a3 = 0.0203;

void forwardKinematics(const double *jpos, double *EE)
{
    Eigen::Matrix<double,6,1> q;
    for (int i = 0; i < 6; i++) {
        q[i] = jpos[i];
    }

    Eigen::Matrix<double,3,1> FKE; FKE.setZero();
    FKE << cos(q[0]) * (a2*cos(q[1]) - a3*cos(q[1]+q[2]) + d4*sin(q[1]+q[2])) - d3*sin(q[0]),
           sin(q[0]) * (a2*cos(q[1]) - a3*cos(q[1]+q[2]) + d4*sin(q[1]+q[2])) + d3*cos(q[0]),
           a3*sin(q[1]+q[2])-a2*sin(q[1])+d4*cos(q[1]+q[2]);

    EE[0] = FKE[0];
    EE[1] = FKE[1];
    EE[2] = FKE[2];
}

void differentialKinematics(const double *jpos, double *JacobEE)
{

    Eigen::Matrix<double,6,1> q;
    for (int i = 0; i < 6; i++) {
        q[i] = jpos[i];
    }

    Eigen::Map<Eigen::Matrix<double,3,6>> J5(JacobEE); J5.setZero();
    J5.col(0) << -sin(q[0])*(-a3*cos(q[1] + q[2]) + d4*sin(q[1] + q[2]) + a2*cos(q[1])) - d3*cos(q[0]),
                 cos(q[0])*(-a3*cos(q[1] + q[2]) + d4*sin(q[1] + q[2]) + a2*cos(q[1])) - d3*sin(q[0]),
                                                                       0;
    J5.col(1) << -cos(q[0])*(-d4*cos(q[1] + q[2]) - a3*sin(q[1] + q[2]) + a2*sin(q[1])),
                 -sin(q[0])*(-d4*cos(q[1] + q[2]) - a3*sin(q[1] + q[2]) + a2*sin(q[1])),
                 -d4*sin(q[1] + q[2]) + a3*cos(q[1] + q[2]) - a2*cos(q[1]);
    J5.col(2) << -cos(q[0])*(-d4*cos(q[1] + q[2]) - a3*sin(q[1] + q[2])),
                 -sin(q[0])*(-d4*cos(q[1] + q[2]) - a3*sin(q[1] + q[2])),
                 -d4*sin(q[1] + q[2]) + a3*cos(q[1] + q[2]);

}


#endif
