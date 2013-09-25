/*
  This is the Linear Solver component of the Amanzi code.
  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)
  Konstantin Lipnikov (lipnikov@lanl.gov)

  Uses NKA as a linear solver.  This is effectively equivalent to GMRES with a
  rolling restart, i.e. vectors fall off the end of the space.

  Usage:
*/

#ifndef AMANZI_NKA_OPERATOR_HH_
#define AMANZI_NKA_OPERATOR_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "exceptions.hh"
#include "VerboseObject.hh"

#include "LinearOperatorDefs.hh"
#include "LinearOperator.hh"
#include "NKA_Base.hh"

namespace Amanzi {
namespace AmanziSolvers {

template<class Matrix, class Vector, class VectorSpace>
class LinearOperatorNKA : public LinearOperator<Matrix, Vector, VectorSpace> {
 public:
  LinearOperatorNKA(const Teuchos::RCP<const Matrix>& m, const Teuchos::RCP<const Matrix>& h) :
      LinearOperator<Matrix, Vector, VectorSpace>(m, h) {
    tol_ = 1e-8;
    overflow_tol_ = 3.0e+50;  // mass of the Universe (J.Hopkins)
    max_itrs_ = 100;
    nka_dim_ = 10;
    criteria_ = LIN_SOLVER_RELATIVE_RHS;
    initialized_ = false;
  }

  void Init(Teuchos::ParameterList& plist);

  Teuchos::RCP<Matrix> Clone() const { return Teuchos::rcp(new LinearOperatorNKA(*this)); }

  int ApplyInverse(const Vector& v, Vector& hv) const;

  // mutators
  void set_tolerance(double tol) { tol_ = tol; }
  void set_max_itrs(int max_itrs) { max_itrs_ = max_itrs; }
  void set_criteria(int criteria) { criteria_ = criteria; }
  void add_criteria(int criteria) { criteria_ |= criteria; }

  // accessors
  double residual() { return residual_; }
  int num_itrs() { return num_itrs_; }

 public:
  Teuchos::RCP<VerboseObject> vo_;

 private:
  using LinearOperator<Matrix, Vector, VectorSpace>::m_;
  using LinearOperator<Matrix, Vector, VectorSpace>::h_;
  using LinearOperator<Matrix, Vector, VectorSpace>::name_;

  int max_itrs_, nka_dim_, criteria_;
  double tol_, nka_tol_, overflow_tol_;
  mutable int num_itrs_;
  mutable double residual_;
  mutable bool initialized_;
};



// Apply the inverse, x <-- A^-1 b
template<class Matrix, class Vector, class VectorSpace>
int LinearOperatorNKA<Matrix, Vector, VectorSpace>::ApplyInverse(const Vector& f, Vector& x) const 
{
  NKA_Base<Vector, VectorSpace> nka(nka_dim_, nka_tol_, f.Map());
  nka.Restart();

  residual_ = 0.0;
  num_itrs_ = 0;

  Teuchos::RCP<Vector> dx  = Teuchos::rcp(new Vector(x));
  Teuchos::RCP<Vector> dxp = Teuchos::rcp(new Vector(x));  // preconditioned correction
  Teuchos::RCP<Vector> r   = Teuchos::rcp(new Vector(x));

  double fnorm, xnorm;
  f.Norm2(&fnorm);
  if (fnorm == 0.0) {
    x.PutScalar(0.0);
    return criteria_;  // Zero solution satifies all criteria.
  }

  x.Norm2(&xnorm);

  m_->Apply(x, *r);  // r = f - A * x
  r->Update(1.0, f, -1.0);

  double rnorm0;
  r->Norm2(&rnorm0);
  residual_ = rnorm0;

  if (initialized_) {
    if (vo_->getVerbLevel() >= Teuchos::VERB_EXTREME) {
      Teuchos::OSTab tab = vo_->getOSTab();
      *(vo_->os()) << num_itrs_ << " ||r||=" << residual_ << endl;
    }
  }
  if (criteria_ == LIN_SOLVER_RELATIVE_RHS) {
    if (rnorm0 < tol_ * fnorm) return LIN_SOLVER_RELATIVE_RHS;
  }

  if (residual_ > overflow_tol_) {
    if (vo_->os_OK(Teuchos::VERB_MEDIUM))
      *vo_->os() << "Overflow: (" << num_itrs_ << " itrs) residual = "
                 << rnorm0 << " (tol=" << overflow_tol_ << ")" << std::endl;
    return LIN_SOLVER_RESIDUAL_OVERFLOW;
  }

  bool done = false;
  while (!done) {
    h_->ApplyInverse(*r, *dxp);
    nka.Correction(*dxp, *dx);
    x.Update(1.0, *dx, 1.0);

    m_->Apply(x, *r);  // r = f - A * x
    r->Update(1.0, f, -1.0);

    double rnorm;
    r->Norm2(&rnorm);
    residual_ = rnorm;

    num_itrs_++;

    if (initialized_) {
      if (vo_->getVerbLevel() >= Teuchos::VERB_EXTREME) {
        Teuchos::OSTab tab = vo_->getOSTab();
        *(vo_->os()) << num_itrs_ << " ||r||=" << residual_ << endl;
      }
    }
    if (rnorm > overflow_tol_) return LIN_SOLVER_RESIDUAL_OVERFLOW;

    // Return the first criterion which is fulfilled.
    if (criteria_ & LIN_SOLVER_RELATIVE_RHS) {
      if (rnorm < tol_ * fnorm) return LIN_SOLVER_RELATIVE_RHS;
    } else if (criteria_ & LIN_SOLVER_RELATIVE_RESIDUAL) {
      if (rnorm < tol_ * rnorm0) return LIN_SOLVER_RELATIVE_RESIDUAL;
    } else if (criteria_ & LIN_SOLVER_ABSOLUTE_RESIDUAL) {
      if (rnorm < tol_) return LIN_SOLVER_ABSOLUTE_RESIDUAL;
    }

    done = num_itrs_ > max_itrs_;
  }

  if (initialized_) {
    if (vo_->os_OK(Teuchos::VERB_MEDIUM))
      *vo_->os() << "Failed (" << num_itrs_ << " itrs) residual = "
                 << residual_ << " (tol=" << tol_ << ")" << std::endl;
  }

  return LIN_SOLVER_MAX_ITERATIONS;
}


/* ******************************************************************
 * Initialization from a parameter list. Available parameters:
 * "error tolerance" [double] default = 1e-6
 * "maximum number of iterations" [int] default = 100
 * "convergence criteria" Array(string) default = "{relative rhs}"
 ****************************************************************** */
template<class Matrix, class Vector, class VectorSpace>
void LinearOperatorNKA<Matrix, Vector, VectorSpace>::Init(Teuchos::ParameterList& plist)
{
  vo_ = Teuchos::rcp(new VerboseObject("Amanzi::NKA_Solver", plist));

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

  // parameters for NKA
  nka_dim_ = plist.get<int>("max nka vectors", 10);
  nka_dim_ = std::min<int>(nka_dim_, max_itrs_);
  nka_tol_ = plist.get<double>("nka vector tolerance", 0.05);
}

} // namespace
} // namespace

#endif
