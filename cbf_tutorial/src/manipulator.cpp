#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <cbf.h>
#include "matplotlibcpp.h"
#include "kinematicsPUMA560.h"

using namespace std;
namespace plt = matplotlibcpp;

const uint32_t nx = 6;
const uint32_t nu = 6;
const uint32_t npSS = 1;
const double lb[nu] = {-1.0,-1.0,-1.0,-1.0,-1.0,-1.0};
const double ub[nu] = {1.0,1.0,1.0,1.0,1.0,1.0};

//double obstacle_pos[3] = {1,1.6,2};
double obstacle_pos[3] = {0,.55,.55};
//double radius = 0.3;
double radius = 0.1;
double wrist_radius = 0.1;

CBF::CBF *cbf;


void safetySet(const double *x, double *h, double *Dh)
{
    Eigen::Map<Eigen::Matrix<double,npSS,1>> hEigen(h);
    Eigen::Map<Eigen::Matrix<double,npSS,nx>> DhEigen(Dh);
    
    Eigen::Matrix<double,nx,1> q;
    for (int i = 0; i < nx; i++) {
        q[i] = x[i];
    }

    double EE[3];
    Eigen::Map<Eigen::Vector3d> E(EE);
    
    forwardKinematics(x, EE);
    
    double JacobEE[3*nx];    
    Eigen::Map<Eigen::Matrix<double,3,nx>> Jpos(JacobEE);

    differentialKinematics(x, JacobEE);
    
    hEigen[0] = (E[0]-obstacle_pos[0])*(E[0]-obstacle_pos[0])+(E[1]-obstacle_pos[1])*(E[1]-obstacle_pos[1])+(E[2]-obstacle_pos[2])*(E[2]-obstacle_pos[2]) - (radius+wrist_radius)*(radius+wrist_radius);

    Eigen::Matrix<double,1,3> dhdx;
    dhdx[0] = 2*(E[0]-obstacle_pos[0]);
    dhdx[1] = 2*(E[1]-obstacle_pos[1]);
    dhdx[2] = 2*(E[2]-obstacle_pos[2]);
    DhEigen = dhdx*Jpos;
}

void dynamics(const double *x, double *f, double *g)
{

    for(uint32_t i=0;i<nx;i++)
    {
        f[i] = 0.0;
    }

    for(uint32_t i=0;i<nx;i++)
    {
        for(uint32_t j=0;j<nu;j++)
        {
            if (i == j)
                g[i+(j*nx)] = 1;
            else
                g[i+(j*nx)] = 0;
        }
    }
}

int main()
{
    CBF::CBF::Options opts;
    opts.relaxLb = 5;

    cbf = new CBF::CBF(nx,nu,npSS,safetySet,dynamics);
    cbf->initialize(lb,ub,opts);

    const double tEnd = 3.0;
    const double dt = 0.001;
    double xNow[6] = {0.0,0.0,0.0,0.0,0.0,0.0};
    double uDesNow[6] = {1.0,0.0,0.0,0.0,0.0,0.0};
    double uActNow[6] = {0.0,0.0,0.0,0.0,0.0,0.0};
    double tNow = 0.0;
    double relax = 0.0;
    int32_t rc;

	// to log the CBF values
	double hNow[npSS];
	double DhNow[npSS*nx];
    
    double EENow[3] = {0.0, 0.0, 0.0};

    int plotLength = (int) (tEnd/dt);
    std::vector<double> hPlot(plotLength);
    std::vector<double> q0Plot(plotLength);
    std::vector<double> q1Plot(plotLength);
    std::vector<double> q2Plot(plotLength);
    std::vector<double> q3Plot(plotLength);
    std::vector<double> q4Plot(plotLength);
    std::vector<double> q5Plot(plotLength);
    std::vector<double> tPlot(plotLength);
    std::vector<double> u0DesPlot(plotLength);
    std::vector<double> u0ActPlot(plotLength);
    std::vector<double> u1DesPlot(plotLength);
    std::vector<double> u1ActPlot(plotLength);
    std::vector<double> u2DesPlot(plotLength);
    std::vector<double> u2ActPlot(plotLength);
    std::vector<double> u3DesPlot(plotLength);
    std::vector<double> u3ActPlot(plotLength);
    std::vector<double> u4DesPlot(plotLength);
    std::vector<double> u4ActPlot(plotLength);
    std::vector<double> u5DesPlot(plotLength);
    std::vector<double> u5ActPlot(plotLength);
    std::vector<double> xEEPlot(plotLength);
    std::vector<double> yEEPlot(plotLength);
    std::vector<double> zEEPlot(plotLength);
    
    int plotCtr = 0;

    bool optionsUpdated = false;
    while(tNow<tEnd)
    {
        if(!optionsUpdated && tNow>(tEnd-2*dt))
        {
            optionsUpdated = true;
            //cbf->updateOptions(opts);
        }

        rc = cbf->filter(xNow,uDesNow,uActNow,relax);
        safetySet(xNow,hNow,DhNow);

        if(rc!=1)
        {
            cout << "CBF failed" << endl;
        }

        // Integrate
        double fCl[nx] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        double f[nx];
        double g[nx*nu];
        dynamics(xNow, f, g);
        for(uint32_t i = 0; i<nx; i++)
        {
            fCl[i]+=f[i];
            for(uint32_t k = 0; k<nu; k++)
            {
                fCl[i]+= uActNow[k]*g[i+k*nx];
            }
        }

        for(uint32_t i = 0; i<nx; i++)
        {
            xNow[i] += dt*fCl[i];
        }
        tNow += dt;

       
        forwardKinematics(xNow, EENow);
        xEEPlot[plotCtr] = EENow[0];
        yEEPlot[plotCtr] = EENow[1];
        zEEPlot[plotCtr] = EENow[2];
        
        hPlot[plotCtr] = hNow[0];
        q0Plot[plotCtr] = xNow[0];
        q1Plot[plotCtr] = xNow[1];
        q2Plot[plotCtr] = xNow[2];
        q3Plot[plotCtr] = xNow[3];
        q4Plot[plotCtr] = xNow[4];
        q5Plot[plotCtr] = xNow[5];
        u0ActPlot[plotCtr] = uActNow[0];
        u0DesPlot[plotCtr] = uDesNow[0];
        u1ActPlot[plotCtr] = uActNow[1];
        u1DesPlot[plotCtr] = uDesNow[1];
        u2ActPlot[plotCtr] = uActNow[2];
        u2DesPlot[plotCtr] = uDesNow[2];
        u3ActPlot[plotCtr] = uActNow[3];
        u3DesPlot[plotCtr] = uDesNow[3];
        u4ActPlot[plotCtr] = uActNow[4];
        u4DesPlot[plotCtr] = uDesNow[4];
        u5ActPlot[plotCtr] = uActNow[5];
        u5DesPlot[plotCtr] = uDesNow[5];
        tPlot[plotCtr] = tNow;
        plotCtr++;
        // Print data
        uint32_t printAcc = 3U;
        cout << std::fixed;
        cout << std::setprecision(3);
        cout << "t: " << tNow ;
        cout << std::scientific;
        cout << std::setprecision(printAcc);
        cout << " | q0:";
        cout << std::setw(printAcc+7) << xNow[0];
        cout << " | q1:";
        cout << std::setw(printAcc+7) << xNow[1];
        cout << " | uDes0:";
        cout << std::setw(printAcc+7) << uDesNow[0];
        cout << " | uAct0:";
        cout << std::setw(printAcc+7) << uActNow[0];
        cout << " | dt(us):";
        cout << std::fixed;
        cout << endl;
    }

    cout << "Finished" << endl;

    plt::backend("Agg");
    plt::named_plot("q0",tPlot,q0Plot);
    plt::named_plot("q1",tPlot,q1Plot);
    plt::named_plot("q2",tPlot,q2Plot);
    plt::xlabel("t (s)");
    plt::ylabel("q (rad)");
    plt::legend();
       cout << "Major joints positions (close figure to proceed)" << endl;
    //plt::show();
    plt::save("/results/manipulator_joints_pos.png");

    plt::clf();
    plt::named_plot("u0Act",tPlot,u0ActPlot,"r");
    plt::named_plot("u0Des",tPlot,u0DesPlot,"--r");
    plt::named_plot("u1Act",tPlot,u1ActPlot,"g");
    plt::named_plot("u1Des",tPlot,u1DesPlot,"--g");
    plt::named_plot("u2Act",tPlot,u2ActPlot,"b");
    plt::named_plot("u2Des",tPlot,u2DesPlot,"--b");
    plt::xlabel("t (s)");
    plt::ylabel("Inputs (rad/s)");
    plt::legend();
       cout << "Major joints velocities (close figure to proceed)" << endl;
    //plt::show();
    plt::save("/results/manipulator_joints_vel.png");
     
    plt::clf(); 
    std::map<std::string, std::string> keywords;
    keywords.insert(std::pair<std::string, std::string>("label", "EE Trajectory") );
    plt::plot3(xEEPlot,yEEPlot,zEEPlot,keywords);
    // Plot spherical obstacle
    std::vector<double> xObs, yObs, zObs;
    for (double u = 0; u <= 2.0*M_PI;  u += M_PI/100) {
        for (double v = 0; v <= 2.0*M_PI;  v += M_PI/100) {
            xObs.push_back(radius * sin(u) * cos(v) + obstacle_pos[0]);
            yObs.push_back(radius * sin(u) * sin(v) + obstacle_pos[1]);
            zObs.push_back(radius * cos(u) + obstacle_pos[2]);
        }
    }
    keywords.clear();
    keywords.insert(std::pair<std::string, std::string>("label", "Obstacle") );
    plt::plot3(xObs,yObs,zObs,keywords);
    plt::xlabel("X");
    plt::ylabel("Y");
    plt::set_zlabel("Z");
    plt::legend();
        cout << "EE motion in 3D (close figure to proceed)" << endl;
    //plt::show();
    plt::save("/results/manipulator_3D_motion.png");

    plt::clf();
    plt::named_plot("h",tPlot,hPlot);
    plt::xlabel("t (s)");
    plt::ylabel("CBF Value");
    plt::legend();
        cout << "CBF Value (close figure to proceed)" << endl;
    //plt::show();
    plt::save("/results/manipulator_cbf.png");
    
    std::ofstream logfile("/results/manipulator_log.csv");
	logfile << "t;q0;q1;q2;q3;q4;q5;u0Act;u1Act;u2Act;u3Act;u4Act;u5Act;u0Des;u1Des;u2Des;u3Des;u4Des;u5Des;xEE;yEE;zEE;h" << endl;
    for(int i=0; i<plotLength; i++)
    {
		logfile << tPlot[i] << ";" << q0Plot[i] << ";" << q1Plot[i] << ";" << q2Plot[i] << ";" 
               << q3Plot[i] << ";" << q4Plot[i] << ";" << q5Plot[i] << ";" 
               << u0ActPlot[i] << ";" << u1ActPlot[i] << ";" << u2ActPlot[i] << ";" 
               << u3ActPlot[i] << ";" << u4ActPlot[i] << ";" << u5ActPlot[i] << ";"
		       << u0DesPlot[i] << ";" << u1DesPlot[i] << ";" << u2DesPlot[i] << ";"
		       << u3DesPlot[i] << ";" << u4DesPlot[i] << ";" << u5DesPlot[i] << ";"
			   << xEEPlot[i] << ";" << yEEPlot[i] << ";" << zEEPlot[i] << ";"
               << hPlot[i] << endl;
    }
    logfile.close();

    return 0;
}
