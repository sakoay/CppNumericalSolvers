// CppNumericalSolver
#ifndef ISOLVER_H_
#define ISOLVER_H_

#include <functional>
#include "isolver.h"
#include "../meta.h"
#include "../problem.h"

namespace cppoptlib {

template<typename T, int Ord>
class ISolver {
public:
    typedef Criteria<T> Criteria;

protected:
    const int order_ = Ord;
    Criteria m_stop, m_current;
    Status m_status = Status::NotStarted;
    DebugLevel m_debug = DebugLevel::None;

public:
    ISolver() {
        m_stop = Criteria::defaults();
        m_current.reset();
    }

    ISolver(const Criteria &s) {
        m_stop = s;
        m_current.reset();
    }

    void setStopCriteria(const Criteria &s) { m_stop = s; }
    const Criteria &criteria() { return m_current; }
    const Status &status() { return m_status; }
    void setDebug(const DebugLevel &d) { m_debug = d; }

    /**
     * @brief minimize an objective function given a gradient (and optinal a hessian)
     * @details this is just the abstract interface
     *
     * @param x0 starting point
     * @param funObjective objective function
     * @param funGradient gradient function
     * @param funcHession hessian function
     */
    virtual void minimize(Problem<T> &objFunc, Vector<T> & x0) = 0;

};

} /* namespace cppoptlib */

#endif /* ISOLVER_H_ */
