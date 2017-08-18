// CppNumericalSolver
#ifndef BFGSSOLVER_H_
#define BFGSSOLVER_H_

#include <iostream>
#include <Eigen/LU>
#include "isolver.h"
#include "../linesearch/morethuente.h"

namespace cppoptlib {

template<typename ProblemType>
class BfgsSolver : public ISolver<ProblemType, 1> {
  public:
    using Superclass = ISolver<ProblemType, 1>;
    using typename Superclass::Scalar;
    using typename Superclass::TVector;
    using typename Superclass::THessian;

    void minimize(ProblemType &objFunc, TVector & x0) {
        const size_t DIM = x0.rows();
        THessian H = THessian::Identity(DIM, DIM);
        TVector grad(DIM);
        TVector x_old = x0;
        this->m_current.reset();
#ifdef SAK_DEBUG_PRINT
        Scalar f_old = objFunc.valueAndGradient(x0, grad);
#else
        objFunc.gradient(x0, grad);
#endif
        if (this->m_stop.xDelta <= 0)
          this->m_stop.xDelta = 1e-7;      // SAK tunable criteria


        do {
            TVector searchDir = -1 * H * grad;
            // check "positive definite"
            Scalar phi = grad.dot(searchDir);

            // positive definit ?
            if (phi > 0 || phi != phi) {
                // no, we reset the hessian approximation
                H = THessian::Identity(DIM, DIM);
                searchDir = -1 * grad;
            }

            const Scalar rate = MoreThuente<ProblemType, 1>::linesearch(x0, searchDir, objFunc) ;
            x0 = x0 + rate * searchDir;

            TVector grad_old = grad;
#ifdef SAK_DEBUG_PRINT
            Scalar f = objFunc.valueAndGradient(x0, grad);
#else
            objFunc.gradient(x0, grad);
#endif
            TVector s = rate * searchDir;
            TVector y = grad - grad_old;

            const Scalar rho = 1.0 / y.dot(s);
            H = H - rho * (s * (y.transpose() * H) + (H * y) * s.transpose()) + rho * rho * (y.dot(H * y) + 1.0 / rho)
                * (s * s.transpose());
            // std::cout << "iter: "<<iter<< " f = " <<  objFunc.value(x0) << " ||g||_inf "<<gradNorm   << std::endl;
#ifdef SAK_DEBUG_PRINT
            mexPrintf(" [%5d]  f = %10.5g (%10.5g), ||g|| = %10.5g, x = [", this->m_current.iterations, f, f_old - f, grad.norm());
            for (Eigen::Index iX = 0; iX < x0.size(); ++iX)
              mexPrintf("%s %8.3g", iX ? "," : "", x0[iX]);
            mexPrintf(" ]\n");
            mexEvalString("drawnow");
            f_old = f;
#endif

            this->m_current.xDelta = (x_old-x0).template lpNorm<Eigen::Infinity>();
            x_old = x0;
            ++this->m_current.iterations;
            this->m_current.gradNorm = grad.template lpNorm<Eigen::Infinity>();
            this->m_status = checkConvergence(this->m_stop, this->m_current);
        } while (objFunc.callback(this->m_current, x0) && (this->m_status == Status::Continue));

    }

};

}
/* namespace cppoptlib */

#endif /* BFGSSOLVER_H_ */
