#include "cbf.h"

namespace CBF
{
    CBF::CBF(const uint32_t nx,
             const uint32_t nu,
             const uint32_t npSS,
             std::function<void(const double* /*x*/,
                                double* /*h*/,
                                double* /*Dh*/)> safetySet,
             std::function<void(const double* /*x*/,
                                double* /*f*/,
                                double* /*g*/)> dynamics,
             const uint32_t npSSmax,
             const QPSOLVER qpSolverType,
             const bool diagonalCost):
            nx_(nx),
            nu_(nu),
            nv_(nu+1),
            npSS_(npSS),
            npSSmax_(npSSmax>npSS?npSS:npSSmax),
            nc_(npSSmax_),
            safetySet_(safetySet),
            dynamics_(dynamics),
            options_()
    {
        H_ = new double[nv_*nv_]{0.0};
        c_ = new double[nv_]{0.0};
        A_ = new double[nc_*nv_]{0.0};
        b_ = new double[nc_]{0.0};
        lb_ = new double[nv_]{0.0};
        ub_ = new double[nv_]{0.0};

        switch(qpSolverType)
        {
            case QPSOLVER::OSQP:
            {
                QPsolver_ = new QPWrapperOsqp(nv_,nc_,diagonalCost);
                break;
            }
            default :
            {
                QPsolver_ = nullptr;
                break;
            }
        }
    }

    CBF::~CBF(void)
    {
        delete[] H_;
        delete[] c_;
        delete[] A_;
        delete[] b_;
        delete[] lb_;
        delete[] ub_;

        if(QPsolver_!=nullptr)
            delete QPsolver_;
    }


    int32_t CBF::initialize(const double lb[],
                            const double ub[])
    {
#ifdef CBF_DEBUG
        std::cout << "Initializing CBF" << std::endl;
#endif

        // Initialize cost hessian
		double *Hidentity;
		Hidentity = new	double[nu_*nu_];
        for(uint32_t j=0; j<nu_;j++)
        {
            for(uint32_t i=0; i<nu_;i++)
            {
                if(i==j)
                    Hidentity[i+(j*nu_)] = 1.0;
                else
                    Hidentity[i+(j*nu_)] = 0.0;
            }
        }
        updateH(Hidentity);
        H_[(nv_-1)+((nv_-1)*nv_)] = options_.relaxCost;

        // Initialize bounds
        memcpy(lb_,lb,nu_*sizeof(double));
        lb_[nu_] = options_.relaxLb;
        memcpy(ub_,ub,nu_*sizeof(double));
        ub_[nu_] = options_.inf;

        // Initialize constraints
		double *x0;
		x0 = new double[nx_];
        for (int i = 0; i < nx_; i++)
            x0[i] = 0;
        updateConstraints(x0);

        // Initialize linear cost
        c_[nu_] = -2.0*options_.relaxCost*options_.relaxLb;
		double *uDes0;
		uDes0 = new double[nu_];
        for (int i = 0; i < nu_; i++)
            uDes0[i] = 0;
        updateCost(uDes0);

        QPsolver_->initialize(H_, c_, A_, b_, lb_, ub_);
        QPsolver_->solve();

#ifdef CBF_DEBUG
        double *solutionFull;
        solutionFull = new double[nv_];
		QPsolver_->getSolution(solutionFull);
#endif

        return 1;
    }

    int32_t CBF::initialize(const double lb[],
                            const double ub[],
                            const Options &options)
    {
        options_ = options;
        return initialize(lb,ub);
    }

    int32_t CBF::filter(const double x[],
                        const double uDes[],
                        double uAct[])
    {
        double relax;
        updateCost(uDes);
        return filter(x,nullptr,c_,uAct,relax);
    }

    int32_t CBF::filter(const double x[],
                        const double uDes[],
                        double uAct[],
                        double &relax)
    {
        updateCost(uDes);
        return filter(x,nullptr,c_,uAct,relax);
    }

    int32_t CBF::filter(const double x[],
                        const double H[],
                        const double c[],
                        double uAct[])
    {
        double relax;
        return filter(x,H,c,uAct,relax);
    }

    int32_t CBF::filter(const double x[],
                        const double H[],
                        const double c[],
                        double uAct[],
                        double &relax)
    {
#ifdef CBF_DEBUG
        std::cout << "Filtering"<< std::endl;
#endif

        updateConstraints(x);
        if(H==nullptr)
            QPsolver_->updateCost(nullptr,c);
        else
        {
            updateH(H);
            QPsolver_->updateCost(H_,c);
        }

        QPsolver_->updateA(A_);
        QPsolver_->updateb(b_);
        int32_t rtCode = QPsolver_->solve();

        if(static_cast<QPWrapperAbstract::SOLVER_STATUS>(rtCode)==QPWrapperAbstract::SOLVER_STATUS::FEASIBLE)
        {
			double *solutionFull;
			solutionFull = new double[nv_];
            QPsolver_->getSolution(solutionFull);
            memcpy(uAct,solutionFull,nu_*sizeof(double));
            inputSaturate(uAct);
            relax = solutionFull[nu_];
            return 1;
        }
        else
            return -1;
    }


    int32_t CBF::updateOptions(const Options &options)
    {
#ifdef CBF_DEBUG
        std::cout << "Updating options"<< std::endl;
#endif

        options_ = options;
        return updateOptions();
    }

    int32_t CBF::updateOptions(void)
    {
        H_[(nv_-1)+((nv_-1)*nv_)] = options_.relaxCost;
        c_[nu_] = -2.0*options_.relaxCost*options_.relaxLb;
        lb_[nu_] = options_.relaxLb;
        QPsolver_->updateBounds(lb_,nullptr);
        QPsolver_->updateCost(H_,c_);
        return 1;
    }

    int32_t CBF::updateConstraints(const double x[])
    {
#ifdef CBF_DEBUG
        std::cout << "Updating constraints"<< std::endl;
#endif

		double *hFull;
		hFull = new double[npSS_];
		double *DhFull;
		DhFull = new double[npSS_*nx_];
		double *f;
		f = new double[nx_];
		double *g;
		g = new double[nx_*nu_];

        safetySet_(x,hFull,DhFull);
        dynamics_(x,f,g);

        double *h;
        double *Dh;

        if(npSSmax_<npSS_)
        {
            h = new double[npSSmax_];
            Dh = new double[npSSmax_*nx_];

            std::vector<uint32_t> idxOfMinH(npSS_);
            for(uint32_t i=0; i<npSS_;i++)
                idxOfMinH[i] = i;

            sort(idxOfMinH.begin(), idxOfMinH.end(), [&hFull](int a, int b) -> bool {return hFull[a]<hFull[b];});
            for(uint32_t i=0; i<npSSmax_; i++)
            {
                h[i] = hFull[idxOfMinH[i]];
                for(uint32_t j=0; j<nx_; j++)
                {
                    Dh[i+j*npSSmax_] = DhFull[idxOfMinH[i]+j*npSS_];
                }
            }
        }
        else
        {
            h = hFull;
            Dh = DhFull;
        }

		double *Lfh;
		Lfh = new double[npSSmax_];
		double *Lgh;
		Lgh = new double[npSSmax_*nu_];

        // Compute Lfh
        matrixVectorMultiply(Dh,npSSmax_,nx_,
                             f,nx_,
                             Lfh);

        // Compute Lgh
        matrixMultiply(Dh,npSSmax_,nx_,
                       g,nx_,nu_,
                       Lgh);

        // Fillup A and b
        for(uint32_t i=0;i<npSSmax_;i++)
        {
            for(uint32_t j=0;j<nu_;j++)
            {
                A_[i+(j*nc_)] = Lgh[i+(j*npSSmax_)];
            }
            A_[i+(nu_*nc_)] = h[i];
            b_[i] = -Lfh[i];
        }

        if(npSSmax_<npSS_)
        {
            delete[] h;
            delete[] Dh;
        }

        return 1;
    }

    int32_t CBF::updateCost(const double uDes[])
    {
#ifdef CBF_DEBUG
        std::cout << "Updating cost"<< std::endl;
#endif

        for(uint32_t i=0; i<nu_;i++)
        {
            c_[i] = -2.0*uDes[i];
        }
        return 1;
    }

    int32_t CBF::updateH(const double H[])
    {
#ifdef CBF_DEBUG
        std::cout << "Building H"<< std::endl;
#endif

        for(uint32_t j=0; j<nu_;j++)
        {
            for(uint32_t i=0; i<nu_;i++)
            {
                H_[i+(j*nv_)] = H[i+(j*nu_)];
            }
        }
        return 1;
    }

    void CBF::inputSaturate(double u[])
    {
        for(uint32_t i = 0; i<nu_; i++)
        {
            if(u[i]>ub_[i])
                u[i] = ub_[i];
            else if(u[i]<lb_[i])
                u[i] = lb_[i];
        }
    }
} //end CBF namespace
