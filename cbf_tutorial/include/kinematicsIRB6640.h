#ifndef _KINEMATICSIRB6640_H_
#define _KINEMATICSIRB6640_H_

#include <functional>
#include <Eigen/Core>

#ifndef M_PI
   #define M_PI 3.14159265358979323846
#endif

// Constants for ABB IRB6640, assuming the following DH parameters:
/*
+-------+--------+-------+-------+-------+
| joint | d      | a     | alpha | offset|
+-------+--------+-------+-------+-------+
|    q0 | 0.78   | 0.32  | -PI/2 | 0     |
|    q1 | 0      | 1.075 | 0     | -PI/2 |
|    q2 | 0      | 0.2   | -PI/2 | 0     |
|    q3 | 1.3925 | 0     | -PI/2 | PI    |
|    q4 | 0      | 0     |  PI/2 | 0     |
|    q5 | 0.2    | 0     | 0     | 0     |
+-------+--------+-------+-------+-------+
*/
const double d0 = 0.780;
const double a0 = 0.320;
const double a1 = 1.075;
const double a2 = 0.2;
const double d3 = 1.392;
const double d5 = 0.2;

void forwardKinematics(const double *jpos, double *EE)
{
    Eigen::Matrix<double,6,1> q;
    for (int i = 0; i < 6; i++) {
        q[i] = jpos[i];
    }

    Eigen::Matrix<double,3,1> FK5; FK5.setZero();
    FK5 << a0*cos(q(0)) + a1*cos(q(0))*sin(q(1)) + a2*(cos(q(0))*cos(q(2))*sin(q(1)) + cos(q(0))*cos(q(1))*sin(q(2))) + d3*(cos(q(0))*cos(q(1))*cos(q(2)) - cos(q(0))*sin(q(1))*sin(q(2))) +
           d5*(cos(q(4))*(cos(q(0))*cos(q(1))*cos(q(2)) - cos(q(0))*sin(q(1))*sin(q(2))) - (cos(q(3))*(cos(q(0))*cos(q(2))*sin(q(1)) + cos(q(0))*cos(q(1))*sin(q(2))) + sin(q(0))*sin(q(3)))*sin(q(4))),
            a0*sin(q(0)) + a1*sin(q(0))*sin(q(1)) + a2*(cos(q(2))*sin(q(0))*sin(q(1)) + cos(q(1))*sin(q(0))*sin(q(2))) + d3*(cos(q(1))*cos(q(2))*sin(q(0)) - sin(q(0))*sin(q(1))*sin(q(2))) +
            d5*(cos(q(4))*(cos(q(1))*cos(q(2))*sin(q(0)) - sin(q(0))*sin(q(1))*sin(q(2))) - (cos(q(3))*(cos(q(2))*sin(q(0))*sin(q(1)) + cos(q(1))*sin(q(0))*sin(q(2))) - cos(q(0))*sin(q(3)))*sin(q(4))),
            d0 + a1*cos(q(1)) + d3*(-(cos(q(2))*sin(q(1))) - cos(q(1))*sin(q(2))) + a2*(cos(q(1))*cos(q(2)) - sin(q(1))*sin(q(2))) +
            d5*(cos(q(4))*(-(cos(q(2))*sin(q(1))) - cos(q(1))*sin(q(2))) - cos(q(3))*(cos(q(1))*cos(q(2)) - sin(q(1))*sin(q(2)))*sin(q(4)));
 
    EE[0] = FK5[0];
    EE[1] = FK5[1];
    EE[2] = FK5[2];
}

void differentialKinematics(const double *jpos, double *JacobEE)
{

    Eigen::Matrix<double,6,1> q;
    for (int i = 0; i < 6; i++) {
        q[i] = jpos[i];
    }

    Eigen::Map<Eigen::Matrix<double,3,6>> J5(JacobEE); J5.setZero();
    J5.col(0) << -(a0*sin(q(0))) - a1*sin(q(0))*sin(q(1)) + a2*(-(cos(q(2))*sin(q(0))*sin(q(1))) - cos(q(1))*sin(q(0))*sin(q(2))) + d3*(-(cos(q(1))*cos(q(2))*sin(q(0))) + sin(q(0))*sin(q(1))*sin(q(2))) +
                 d5*(cos(q(4))*(-(cos(q(1))*cos(q(2))*sin(q(0))) + sin(q(0))*sin(q(1))*sin(q(2))) -
                     (cos(q(3))*(-(cos(q(2))*sin(q(0))*sin(q(1))) - cos(q(1))*sin(q(0))*sin(q(2))) + cos(q(0))*sin(q(3)))*sin(q(4))),
            a0*cos(q(0)) + a1*cos(q(0))*sin(q(1)) + a2*(cos(q(0))*cos(q(2))*sin(q(1)) + cos(q(0))*cos(q(1))*sin(q(2))) + d3*(cos(q(0))*cos(q(1))*cos(q(2)) - cos(q(0))*sin(q(1))*sin(q(2))) +
            d5*(cos(q(4))*(cos(q(0))*cos(q(1))*cos(q(2)) - cos(q(0))*sin(q(1))*sin(q(2))) - (cos(q(3))*(cos(q(0))*cos(q(2))*sin(q(1)) + cos(q(0))*cos(q(1))*sin(q(2))) + sin(q(0))*sin(q(3)))*sin(q(4))),
            0;
    J5.col(1) << a1*cos(q(0))*cos(q(1)) + d3*(-(cos(q(0))*cos(q(2))*sin(q(1))) - cos(q(0))*cos(q(1))*sin(q(2))) + a2*(cos(q(0))*cos(q(1))*cos(q(2)) - cos(q(0))*sin(q(1))*sin(q(2))) +
                 d5*(cos(q(4))*(-(cos(q(0))*cos(q(2))*sin(q(1))) - cos(q(0))*cos(q(1))*sin(q(2))) - cos(q(3))*(cos(q(0))*cos(q(1))*cos(q(2)) - cos(q(0))*sin(q(1))*sin(q(2)))*sin(q(4))),
            a1*cos(q(1))*sin(q(0)) + d3*(-(cos(q(2))*sin(q(0))*sin(q(1))) - cos(q(1))*sin(q(0))*sin(q(2))) + a2*(cos(q(1))*cos(q(2))*sin(q(0)) - sin(q(0))*sin(q(1))*sin(q(2))) +
            d5*(cos(q(4))*(-(cos(q(2))*sin(q(0))*sin(q(1))) - cos(q(1))*sin(q(0))*sin(q(2))) - cos(q(3))*(cos(q(1))*cos(q(2))*sin(q(0)) - sin(q(0))*sin(q(1))*sin(q(2)))*sin(q(4))),
            -(a1*sin(q(1))) + a2*(-(cos(q(2))*sin(q(1))) - cos(q(1))*sin(q(2))) + d3*(-(cos(q(1))*cos(q(2))) + sin(q(1))*sin(q(2))) +
            d5*(cos(q(4))*(-(cos(q(1))*cos(q(2))) + sin(q(1))*sin(q(2))) - cos(q(3))*(-(cos(q(2))*sin(q(1))) - cos(q(1))*sin(q(2)))*sin(q(4)));
    J5.col(2) << d3*(-(cos(q(0))*cos(q(2))*sin(q(1))) - cos(q(0))*cos(q(1))*sin(q(2))) + a2*(cos(q(0))*cos(q(1))*cos(q(2)) - cos(q(0))*sin(q(1))*sin(q(2))) +
                 d5*(cos(q(4))*(-(cos(q(0))*cos(q(2))*sin(q(1))) - cos(q(0))*cos(q(1))*sin(q(2))) - cos(q(3))*(cos(q(0))*cos(q(1))*cos(q(2)) - cos(q(0))*sin(q(1))*sin(q(2)))*sin(q(4))),
            d3*(-(cos(q(2))*sin(q(0))*sin(q(1))) - cos(q(1))*sin(q(0))*sin(q(2))) + a2*(cos(q(1))*cos(q(2))*sin(q(0)) - sin(q(0))*sin(q(1))*sin(q(2))) +
            d5*(cos(q(4))*(-(cos(q(2))*sin(q(0))*sin(q(1))) - cos(q(1))*sin(q(0))*sin(q(2))) - cos(q(3))*(cos(q(1))*cos(q(2))*sin(q(0)) - sin(q(0))*sin(q(1))*sin(q(2)))*sin(q(4))),
            a2*(-(cos(q(2))*sin(q(1))) - cos(q(1))*sin(q(2))) + d3*(-(cos(q(1))*cos(q(2))) + sin(q(1))*sin(q(2))) +
            d5*(cos(q(4))*(-(cos(q(1))*cos(q(2))) + sin(q(1))*sin(q(2))) - cos(q(3))*(-(cos(q(2))*sin(q(1))) - cos(q(1))*sin(q(2)))*sin(q(4)));
    J5.col(3) << -(d5*(cos(q(3))*sin(q(0)) - (cos(q(0))*cos(q(2))*sin(q(1)) + cos(q(0))*cos(q(1))*sin(q(2)))*sin(q(3)))*sin(q(4))),
            -(d5*(-(cos(q(0))*cos(q(3))) - (cos(q(2))*sin(q(0))*sin(q(1)) + cos(q(1))*sin(q(0))*sin(q(2)))*sin(q(3)))*sin(q(4))),d5*(cos(q(1))*cos(q(2)) - sin(q(1))*sin(q(2)))*sin(q(3))*sin(q(4));
    J5.col(4) << d5*(-(cos(q(4))*(cos(q(3))*(cos(q(0))*cos(q(2))*sin(q(1)) + cos(q(0))*cos(q(1))*sin(q(2))) + sin(q(0))*sin(q(3)))) -
                     (cos(q(0))*cos(q(1))*cos(q(2)) - cos(q(0))*sin(q(1))*sin(q(2)))*sin(q(4))),
            d5*(-(cos(q(4))*(cos(q(3))*(cos(q(2))*sin(q(0))*sin(q(1)) + cos(q(1))*sin(q(0))*sin(q(2))) - cos(q(0))*sin(q(3)))) -
                (cos(q(1))*cos(q(2))*sin(q(0)) - sin(q(0))*sin(q(1))*sin(q(2)))*sin(q(4))),
            d5*(-(cos(q(3))*cos(q(4))*(cos(q(1))*cos(q(2)) - sin(q(1))*sin(q(2)))) - (-(cos(q(2))*sin(q(1))) - cos(q(1))*sin(q(2)))*sin(q(4)));

}


#endif
