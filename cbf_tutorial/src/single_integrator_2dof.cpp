#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <cbf.h>
#include "matplotlibcpp.h"

#ifndef M_PI
   #define M_PI 3.14159265358979323846
#endif

using namespace std;
namespace plt = matplotlibcpp;

const uint32_t nx = 2;
const uint32_t nu = 2;
const uint32_t npSS = 1;
const double lb[nu] = { -1.0, -1.0 };
const double ub[nu] = { 1.0, 1.0 };

double obstacle_pos[2] = {1.0, 1.0};
double radius = 0.3;

double target_pos[2] = { 2.0, 2.5 };
double K = 2.5;

CBF::CBF *cbf;

void safetySet(const double *x, double *h, double *Dh)
{
     h[0] = (x[0] - obstacle_pos[0])*(x[0] - obstacle_pos[0]) + (x[1] - obstacle_pos[1])*(x[1] - obstacle_pos[1]) - radius*radius;

     Dh[0] = 2 * (x[0] - obstacle_pos[0]);
     Dh[1] = 2 * (x[1] - obstacle_pos[1]);
     
}

void dynamics(const double *x, double *f, double *g)
{

     const double A[nx*nx] = { 0.0 , 0.0, 0.0, 0.0};
     const double b[nx*nu] = { 1.0 , 0.0, 0.0, 1.0};

     for (uint32_t i = 0;i<nx;i++)
     {
        f[i] = 0.0;
        for (uint32_t j = 0;j<nx;j++)
        {
             f[i] += A[i + (j*nx)] * x[j];
        }
     }

     for (uint32_t i = 0;i<nx;i++)
     {
        for (uint32_t j = 0;j<nu;j++)
        {
             g[i + (j*nx)] = b[i + (j*nx)];
        }
     }
}

int main()
{
     CBF::CBF::Options opts;
     opts.relaxLb = 1.0;

     cbf = new CBF::CBF(nx, nu, npSS, safetySet, dynamics);
     cbf->initialize(lb, ub, opts);

     const double tEnd = 8.0;
     const double dt = 0.001;
     double xNow[nx] = { 0.0, 0.0 };
     double uDesNow[nu] = { 0.0, 0.0 };
     double uActNow[nu] = { 0.0, 0.0 };
     double tNow = 0.0;
     double relax = 0.0;
     int32_t rc;
     
         // to log the CBF values
         double hNow[npSS];
         double DhNow[npSS*nx];

     int plotLength = (int)(tEnd / dt);
     std::vector<double> hPlot(plotLength);
     std::vector<double> x1Plot(plotLength);
     std::vector<double> x2Plot(plotLength);
     std::vector<double> tPlot(plotLength);
     std::vector<double> u1DesPlot(plotLength);
     std::vector<double> u2DesPlot(plotLength);
     std::vector<double> u1ActPlot(plotLength);
     std::vector<double> u2ActPlot(plotLength);
     int plotCtr = 0;

     bool optionsUpdated = false;
     while (tNow<tEnd)
     {
        if (!optionsUpdated && tNow>(tEnd / 2))
        {
             optionsUpdated = true;
             //cbf->updateOptions(opts);
        }

        uDesNow[0] = K*(target_pos[0] - xNow[0]);
        uDesNow[1] = K*(target_pos[1] - xNow[1]);

        rc = cbf->filter(xNow, uDesNow, uActNow, relax);
        safetySet(xNow,hNow,DhNow);

        if (rc != 1)
        {
             cout << "CBF failed" << endl;
        }

        // Integrate
        double fCl[nx] = { 0.0 };
        double f[nx];
        double g[nx*nu];
        dynamics(xNow, f, g);
        for (uint32_t i = 0; i<nx; i++)
        {
             fCl[i] += f[i];
             for (uint32_t k = 0; k<nu; k++)
             {
                fCl[i] += uActNow[k] * g[i + k*nx];
             }
        }

        for (uint32_t i = 0; i<nx; i++)
        {
             xNow[i] += dt*fCl[i];
        }
        tNow += dt;

        hPlot[plotCtr] = hNow[0];
        x1Plot[plotCtr] = xNow[0];
        x2Plot[plotCtr] = xNow[1];
        u1ActPlot[plotCtr] = uActNow[0];
        u2ActPlot[plotCtr] = uActNow[1];
        u1DesPlot[plotCtr] = uDesNow[0];
        u2DesPlot[plotCtr] = uDesNow[1];
        tPlot[plotCtr] = tNow;
        plotCtr++;
        // Print data
        uint32_t printAcc = 3U;
        cout << std::fixed;
        cout << std::setprecision(3);
        cout << "t: " << tNow;
        cout << std::scientific;
        cout << std::setprecision(printAcc);
        cout << " | x1:";
        cout << std::setw(printAcc + 7) << xNow[0];
        cout << " | x2:";
        cout << std::setw(printAcc + 7) << xNow[1];
        cout << " | uDes:";
        cout << std::setw(printAcc + 7) << uDesNow[1];
        cout << " | uAct:";
        cout << std::setw(printAcc + 7) << uActNow[1];
        cout << " | dt(us):";
        cout << std::fixed;
        cout << endl;
     }

     cout << "Finished" << endl;

     plt::backend("Agg");
     plt::plot(x1Plot, x2Plot);
     std::vector<double> dirVectX1(2);
     std::vector<double> dirVectX2(2);
     std::vector<double> quivX1{ x1Plot[0], x1Plot[x1Plot.size() / 3 - 1] };
     std::vector<double> quivX2{ x2Plot[0], x2Plot[x2Plot.size() / 3 - 1] };
     dirVectX1[0] = u1ActPlot[0];
     dirVectX1[1] = u1ActPlot[x1Plot.size() / 3 - 1];
     dirVectX2[0] = u2ActPlot[0];
     dirVectX2[1] = u2ActPlot[x2Plot.size() / 3 - 1];
     plt::quiver(quivX1, quivX2, dirVectX1, dirVectX2);
     // Plot spherical obstacle
     std::vector<double> xObs, yObs;
     for (double u = 0; u <= 2.0*M_PI;  u += M_PI/10) {
          xObs.push_back(radius * cos(u) + obstacle_pos[0]);
          yObs.push_back(radius * sin(u) + obstacle_pos[1]);
     }
     plt::plot(xObs,yObs); 
     plt::xlabel("x1 (m)");
     plt::ylabel("x2 (m)");
     plt::axis("equal");
          cout << "State space plot (close figure to proceed)" << endl;
     //plt::show();
     plt::save("/results/single_integrator_2dof_phase.png");

     plt::clf();
     plt::named_plot("u1Act", tPlot, u1ActPlot);
     plt::named_plot("u1Des", tPlot, u1DesPlot);
     plt::named_plot("u2Act", tPlot, u2ActPlot);
     plt::named_plot("u2Des", tPlot, u2DesPlot);
     plt::xlabel("t (s)");
     plt::ylabel("Inputs (m/s)");
     plt::legend();
         cout << "Actual vs Desired inputs versus time (close figure to proceed)" << endl;
     //plt::show();
     plt::save("/results/single_integrator_2dof_inputs.png");
     

     plt::clf(); 
     plt::named_plot("h",tPlot,hPlot);
     plt::xlabel("t (s)");
     plt::ylabel("CBF Value");
     plt::legend();
     cout << "CBF Value (close figure to proceed)" << endl;
     //plt::show();
     plt::save("/results/single_integrator_2dof_cbf.png");
         
     std::ofstream logfile("/results/single_integrator_2dof_log.csv");
     logfile << "t;x1;x2;u1Act;u2Act;u1Des;u2Des;h" << endl;
     for (int i = 0; i<plotLength; i++)
     {
        logfile << tPlot[i] << ";" << x1Plot[i] << ";" << x2Plot[i] << ";"
             << u1ActPlot[i] << ";" << u2ActPlot[i] << ";"
             << u1DesPlot[i] << ";" << u2DesPlot[i] << ";" 
             << hPlot[i] << endl;
     }
     logfile.close();

     return 0;
}
