#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <cbf.h>
#include "matplotlibcpp.h"

using namespace std;
namespace plt = matplotlibcpp;

const uint32_t nx = 1;
const uint32_t nu = 1;
const uint32_t npSS = 2;
const double lb[nu] = {-1.0};
const double ub[nu] = {1.0};
const double xBound[2] = {-1.0,1.0};
const double vBound[2] = {-1.0,1.0};

CBF::CBF *cbf;

void safetySet(const double *x, double *h, double *Dh)
{
    /*if(x[0]>0)
    {
        h[0] = xBound[1] - x[0];   Dh[0] = -1.0;
    }
    else
    {
        h[0] = x[0] - xBound[0];   Dh[0] = 1.0;
    }*/
	h[0] = xBound[1] - x[0];   Dh[0] = -1.0;
	h[1] = x[0] - xBound[0];   Dh[1] = 1.0;

	/*h[0] = 0.25*(xBound[0] - xBound[1])*(xBound[0] - xBound[1]) - (x[0] - (xBound[0] + xBound[1]) / 2)*(x[0] - (xBound[0] + xBound[1]) / 2);
	Dh[0] = xBound[0] + xBound[1] - 2 * x[0];*/

}

void dynamics(const double *x, double *f, double *g)
{
    const double A[nx*nx] = {0.0};
    const double b[nx*nu] = {1.0};

    for(uint32_t i=0;i<nx;i++)
    {
        f[i] = 0.0;
        for(uint32_t j=0;j<nx;j++)
        {
            f[i] += A[i+(j*nx)]*x[j];
        }
    }

    for(uint32_t i=0;i<nx;i++)
    {
        for(uint32_t j=0;j<nu;j++)
        {
            g[i+(j*nx)] = b[i+(j*nx)];
        }
    }
}

int main()
{
    CBF::CBF::Options opts;
    opts.relaxLb = 2.0;

    cbf = new CBF::CBF(nx,nu,npSS,safetySet,dynamics);
    cbf->initialize(lb,ub,opts);

    const double tEnd = 5.0;
    const double dt = 0.001;
    double xNow[nx] = {0};
    double uDesNow[nu] = {1};
    double uActNow[nu] = {0.0};
    double tNow = 0.0;
    double relax = 0.0;
    int32_t rc;

    int plotLength = (int) (tEnd/dt);
    std::vector<double> xPlot(plotLength);
    std::vector<double> vPlot(plotLength);
    std::vector<double> tPlot(plotLength);
    std::vector<double> uDesPlot(plotLength);
    std::vector<double> uActPlot(plotLength);
    int plotCtr = 0;

    bool optionsUpdated = false;
    while(tNow<tEnd)
    {
        if(!optionsUpdated && tNow>(tEnd/2))
        {
            optionsUpdated = true;
            //cbf->updateOptions(opts);
        }

        rc = cbf->filter(xNow,uDesNow,uActNow,relax);

        if(rc!=1)
        {
            cout << "CBF failed" << endl;
        }

        // Integrate
        double fCl[nx] = {0.0};
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

        xPlot[plotCtr] = xNow[0];
        vPlot[plotCtr] = xNow[1];
        uActPlot[plotCtr] = uActNow[0];
        uDesPlot[plotCtr] = uDesNow[0];
        tPlot[plotCtr] = tNow;
        plotCtr++;
        // Print data
        uint32_t printAcc = 3U;
        cout << std::fixed;
        cout << std::setprecision(3);
        cout << "t: " << tNow ;
        cout << std::scientific;
        cout << std::setprecision(printAcc);
        cout << " | x:";
        cout << std::setw(printAcc+7) << xNow[0];
        cout << " | uDes:";
        cout << std::setw(printAcc+7) << uDesNow[0];
        cout << " | uAct:";
        cout << std::setw(printAcc+7) << uActNow[0];
        cout << " | dt(us):";
        cout << std::fixed;
        cout << endl;
    }

    cout << "Finished" << endl;

    plt::backend("Agg");
    plt::plot(xPlot,uActPlot);
    std::vector<double> dirVectX(2);
    std::vector<double> dirVectV(2);
    std::vector<double> quivX{xPlot[0], xPlot[xPlot.size()/3-1]};
    std::vector<double> quivV{uActPlot[0], uActPlot[uActPlot.size()/3-1]};
    dirVectX[0] = xPlot[1]-xPlot[0];
    dirVectX[1] = xPlot[xPlot.size()/3] - xPlot[xPlot.size()/3-1];
    dirVectV[0] = uActPlot[1]-uActPlot[0];
    dirVectV[1] = uActPlot[uActPlot.size()/3] - uActPlot[uActPlot.size()/3-1];
    plt::quiver(quivX,quivV,dirVectX,dirVectV);
    plt::xlabel("x (m)");
    plt::ylabel("v (m/s)");
        cout << "Phase plot (close figure to proceed)" << endl;
    //plt::show();
    //plt::show();
    plt::save("/results/single_integrator_phase.png");
    
    plt::clf();
    plt::named_plot("uAct",tPlot,uActPlot);
    plt::named_plot("uDes",tPlot,uDesPlot);
    plt::xlabel("t (s)");
    plt::ylabel("Input (m/s)");
    plt::legend();
        cout << "Actual vs Desired input versus time (close figure to proceed)" << endl;
    //plt::show();
    //plt::show();
    plt::save("/results/single_integrator_input.png");
    
    std::ofstream logfile("/results/single_integrator_log.csv");
	logfile << "t;x;v;uAct;uDes" << endl;
    for(int i=0; i<plotLength; i++)
    {
		logfile << tPlot[i] << ";" << xPlot[i] << ";" << vPlot[i] << ";" << uActPlot[i] << ";" << uDesPlot[i] << endl;
    }
    logfile.close();
    
    return 0;
}
