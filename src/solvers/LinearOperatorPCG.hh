/*
  This is the Linear Solver component of the Amanzi code.
  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)
  Konstantin Lipnikov (lipnikov@lanl.gov)

  Conjugate gradient method.
  Usage:
*/

#ifndef AMANZI_PCG_OPERATOR_HH_
#define AMANZI_PCG_OPERATOR_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "VerboseObject.hh"

#include "LinearOperatorDefs.hh"
#include "LinearOperator.hh"

namespace Amanzi {
namespace AmanziSolvers {

template<class Matrix, class Vector, class VectorSpace>
class LinearOperatorPCG : public LinearOperator<Matrix, Vector, VectorSpace> {
 public:
  LinearOperatorPCG(const Teuchos::RCP<const Matrix>& m,
                    const Teuchos::RCP<const Matrix>& h) :
      LinearOperator<Matrix, Vector, VectorSpace>(m, h),
      tol_(1e-6),
      overflow_tol_(3.0e+50),  // mass of the Universe (J.Hopkins)
      max_itrs_(100),
      criteria_(LIN_SOLVER_RELATIVE_RHS),
      initialized_(false) {}

  LinearOperatorPCG(const LinearOperatorPCG& other) :
      LinearOperator<Matrix,Vector,VectorSpace>(other),
      tol_(other.tol_),
      overflow_tol_(other.overflow_tol_),
      max_itrs_(other.max_itrs_),
      num_itrs_(other.num_itrs_),
      residual_(other.residual_),
      criteria_(other.criteria_),
      initialized_(other.initialized_) {}

  virtual Teuchos::RCP<Matrix> Clone() const {
    return Teuchos::rcp(new LinearOperatorPCG(*this)); }

  void Init(Teuchos::ParameterList& plist);

  int ApplyInverse(const Vector& v, Vector& hv) const {
    int ierr = PCG_(v, hv, tol_, max_itrs_, criteria_);
    return ierr;
    // return (ierr > 0) ? 0 : 1;
  }

  // access members
  void set_tolerance(double tol) { tol_ = tol; }
  void set_max_itrs(int max_itrs) { max_itrs_ = max_itrs; }
  void set_criteria(int criteria) { criteria_ = criteria; }
  void add_criteria(int criteria) { criteria_ |= criteria; }
  void set_overflow(double tol) { overflow_tol_ = tol; }

  double residual() { return residual_; }
  int num_itrs() { return num_itrs_; }

 public:
  Teuchos::RCP<VerboseObject> vo_;

 private:
  int PCG_(const Vector& f, Vector& x, double tol, int max_itrs, int criteria) const;

 private:
  using LinearOperator<Matrix, Vector, VectorSpace>::m_;
  using LinearOperator<Matrix, Vector, VectorSpace>::h_;
  using LinearOperator<Matrix, Vector, VectorSpace>::name_;

  int max_itrs_, criteria_;
  double tol_, overflow_tol_;
  mutable int num_itrs_;
  mutable double residual_;
  mutable bool initialized_;
};


/* ******************************************************************
 * PCG input/output data:
 *  f [input]         the right-hand side
 *  x [input/output]  initial guess / final solution
 *  tol [input]       convergence tolerance
 *  max_itrs [input]  maximum number of iterations
 *  criteria [input]  sum of termination critaria
 *
 *  Return value. If it is positive, it indicates the sucessful
 *  convergence criterion (criteria in a few exceptional cases) that
 *  was checked first. If it is negative, it indicates a failure, see
 *  LinearSolverDefs.hh for the error explanation.
 ***************************************************************** */
template<class Matrix, class Vector, class VectorSpace>
int LinearOperatorPCG<Matrix, Vector, VectorSpace>::PCG_(
    const Vector& f, Vector& x, double tol, int max_itrs, int criteria) const
{
  Vector r(f), p(f), v(f);  // construct empty vectors
  num_itrs_ = 0;

  double fnorm;
  f.Norm2(&fnorm);
  if (fnorm == 0.0) {
    x.PutScalar(0.0);
    return criteria;  // Convergence for all criteria
  }

  m_->Apply(x, r);  // r = f - M * x
  r.Update(1.0, f, -1.0);
  double rnorm0;
  r.Norm2(&rnorm0);
  residual_ = rnorm0;

  h_->ApplyInverse(r, p);  // gamma = (H r,r)
  double gamma0;
  p.Dot(r, &gamma0);
  if (gamma0 < 0) return LIN_SOLVER_NON_SPD_APPLY_INVERSE;

  // Ignore all criteria if one iteration is enforced.
  if (rnorm0 > overflow_tol_) return LIN_SOLVER_RESIDUAL_OVERFLOW;

  if (! (criteria & LIN_SOLVER_MAKE_ONE_ITERATION)) {
    if (criteria & LIN_SOLVER_RELATIVE_RHS) {
      if (rnorm0 < tol * fnorm) return LIN_SOLVER_RELATIVE_RHS;
    } else if (criteria & LIN_SOLVER_ABSOLUTE_RESIDUAL) {
      if (rnorm0 < tol) return LIN_SOLVER_ABSOLUTE_RESIDUAL;
    }
  }

  for (int i = 0; i < max_itrs; i++) {
    m_->Apply(p, v);
    double alpha;
    v.Dot(p, &alpha);

    if (alpha < 0.0) return LIN_SOLVER_NON_SPD_APPLY;
    alpha = gamma0 / alpha;

    x.Update( alpha, p, 1.0);
    r.Update(-alpha, v, 1.0);

    h_->ApplyInverse(r, v);  // gamma1 = (H r, r)
    double gamma1;
    v.Dot(r, &gamma1);
    if (gamma1 < 0.0) return LIN_SOLVER_NON_SPD_APPLY_INVERSE;

    double rnorm;
    r.Norm2(&rnorm);
    residual_ = rnorm;

    if (initialized_) {
      if (vo_->getVerbLevel() >= Teuchos::VERB_EXTREME) {
        Teuchos::OSTab tab = vo_->getOSTab();
        *(vo_->os()) << i << " ||r||=" << residual_ << endl;
      }
    }
    if (rnorm > overflow_tol_) return LIN_SOLVER_RESIDUAL_OVERFLOW;

    // Return the first criterion which is fulfilled.
    num_itrs_ = i + 1;
    if (criteria & LIN_SOLVER_RELATIVE_RHS) {
      if (rnorm < tol * fnorm) return LIN_SOLVER_RELATIVE_RHS;
    } else if (criteria & LIN_SOLVER_RELATIVE_RESIDUAL) {
      if (rnorm < tol * rnorm0) return LIN_SOLVER_RELATIVE_RESIDUAL;
    } else if (criteria & LIN_SOLVER_ABSOLUTE_RESIDUAL) {
      if (rnorm < tol) return LIN_SOLVER_ABSOLUTE_RESIDUAL;
    }

    double beta = gamma1 / gamma0;
    gamma0 = gamma1;

    p.Update(1.0, v, beta);
  }

  return LIN_SOLVER_MAX_ITERATIONS;
};


/* ******************************************************************
 * Initialization from a parameter list. Available parameters:
 * "error tolerance" [double] default = 1e-6
 * "maximum number of iterations" [int] default = 100
 * "convergence criteria" Array(string) default = "{relative rhs}"
 ****************************************************************** */
template<class Matrix, class Vector, class VectorSpace>
void LinearOperatorPCG<Matrix, Vector, VectorSpace>::Init(Teuchos::ParameterList& plist)
{
  vo_ = Teuchos::rcp(new VerboseObject("Solvers::PCG", plist));

  tol_ = plist.get<double>("error tolerance", 1e-6);
  max_itrs_ = plist.get<int>("maximum number of iterations", 100);
  overflow_tol_ = plist.get<double>("overflow tolerance", 3.0e+50);

  int criteria(0);
  if (plist.isParameter("convergence criteria")) {
    std::vector<std::string> names;
    names = plist.get<Teuchos::Array<std::string> > ("convergence criteria").toVector();

    for (int i = 0; i < names.size(); i++) {
      if (names[i] == "relative rhs") {
        criteria += LIN_SOLVER_RELATIVE_RHS;
      } else if (names[i] == "relative residual") {
        criteria += LIN_SOLVER_RELATIVE_RESIDUAL;
      } else if (names[i] == "absolute residual") {
        criteria += LIN_SOLVER_ABSOLUTE_RESIDUAL;
      } else if (names[i] == "make one iteration") {
        criteria += LIN_SOLVER_MAKE_ONE_ITERATION;
      }
    }
  } else {
    criteria = LIN_SOLVER_RELATIVE_RHS;
  }

  set_criteria(criteria);
  initialized_ = true;
}

}  // namespace AmanziSolvers
}  // namespace Amanzi

#endif
