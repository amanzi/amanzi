/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Alicia Klinvex (amklinv@sandia.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>

#ifndef AMANZI_BELOS_GMRES_OPERATOR_HH_
#define AMANZI_BELOS_GMRES_OPERATOR_HH_

#include <cmath>

#include "BelosTypes.hpp"
#include "AmanziBelosOPWrapper.hh"
#include "AmanziBelosMVWrapper.hh"
#include "BelosPseudoBlockGmresSolMgr.hpp"

#include "Teuchos_RCP.hpp"
#include "errors.hh"
#include "VerboseObject.hh"

#include "DenseMatrix.hh"
#include "LinearOperator.hh"
#include "LinearOperatorDefs.hh"

namespace Amanzi {
namespace AmanziSolvers {

template <class Matrix, class Vector, class VectorSpace>
class LinearOperatorBelosGMRES
  : public LinearOperator<Matrix, Vector, VectorSpace> {
 public:
  LinearOperatorBelosGMRES(const Teuchos::RCP<const Matrix>& m,
                           const Teuchos::RCP<const Matrix>& h)
    : LinearOperator<Matrix, Vector, VectorSpace>(m, h),
      tol_(1e-6),
      overflow_tol_(3.0e+50), // mass of the Universe (J.Hopkins)
      max_itrs_(100),
      krylov_dim_(10),
      criteria_(LIN_SOLVER_RELATIVE_RHS),
      initialized_(false)
  {}

  void Init(Teuchos::ParameterList& plist);
  void Init() { LinearOperator<Matrix, Vector, VectorSpace>::Init(); }

  int ApplyInverse(const Vector& v, Vector& hv) const;

  // access members
  void set_tolerance(double tol) { tol_ = tol; }
  void set_max_itrs(int max_itrs) { max_itrs_ = max_itrs; }
  void set_criteria(int criteria) { criteria_ = criteria; }
  void add_criteria(int criteria) { criteria_ |= criteria; }
  void set_krylov_dim(int n) { krylov_dim_ = n; }
  void set_overflow(double tol) { overflow_tol_ = tol; }

  double residual() { return residual_; }
  int num_itrs() { return num_itrs_; }
  int returned_code() { return returned_code_; }

 private:
  // using LinearOperator<Matrix, Vector, VectorSpace>::m_; // solving mx=f
  // using LinearOperator<Matrix, Vector, VectorSpace>::h_; // h is
  // preconditioner using LinearOperator<Matrix, Vector, VectorSpace>::name_;

  Teuchos::RCP<VerboseObject> vo_;

  int max_itrs_, criteria_, krylov_dim_;
  double tol_, overflow_tol_;
  mutable int num_itrs_, returned_code_;
  mutable double residual_;
  mutable bool initialized_;
};


/* ******************************************************************
 * Initialization from a parameter list. Available parameters:
 * "error tolerance" [double] default = 1e-6
 * "maximum number of iterations" [int] default = 100
 * "convergence criteria" Array(string) default = "{relative rhs}"
 ****************************************************************** */
template <class Matrix, class Vector, class VectorSpace>
void
LinearOperatorBelosGMRES<Matrix, Vector, VectorSpace>::Init(
  Teuchos::ParameterList& plist)
{
  vo_ = Teuchos::rcp(new VerboseObject("Solvers::BelosGMRES", plist));

  tol_ = plist.get<double>("error tolerance", 1e-16);
  max_itrs_ = plist.get<int>("maximum number of iterations", 100);
  krylov_dim_ = plist.get<int>("size of Krylov space", 10);
  overflow_tol_ = plist.get<double>("overflow tolerance", 3.0e+50);

  int criteria(0);
  if (plist.isParameter("convergence criteria")) {
    std::vector<std::string> names;
    names =
      plist.get<Teuchos::Array<std::string>>("convergence criteria").toVector();

    // TODO: CUSTOM STATUS TESTS
    for (int i = 0; i < names.size(); i++) {
      if (names[i] == "relative rhs") {
        criteria += LIN_SOLVER_RELATIVE_RHS;
      } else if (names[i] == "relative residual") {
        criteria += LIN_SOLVER_RELATIVE_RESIDUAL;
      } else if (names[i] == "absolute residual") {
        criteria += LIN_SOLVER_ABSOLUTE_RESIDUAL;
      } else if (names[i] == "make one iteration") {
        criteria += LIN_SOLVER_MAKE_ONE_ITERATION;
      } else {
        Errors::Message msg;
        msg << "LinearOperatorGMRES: \"convergence criteria\" type \""
            << names[i] << "\" is not recognized.";
        Exceptions::amanzi_throw(msg);
      }
    }
  } else {
    criteria = LIN_SOLVER_RELATIVE_RHS;
  }

  set_criteria(criteria);

  initialized_ = true;
}


template <class Matrix, class Vector, class VectorSpace>
int
LinearOperatorBelosGMRES<Matrix, Vector, VectorSpace>::ApplyInverse(
  const Vector& v, Vector& hv) const
{
  typedef Belos::
    LinearProblem<double, Belos::MultiVec<double>, Belos::Operator<double>>
      LinearProblem;
  typedef Belos::PseudoBlockGmresSolMgr<double,
                                        Belos::MultiVec<double>,
                                        Belos::Operator<double>>
    GmresSolver;

  if (!initialized_) {
    Errors::Message msg("LinearOperatorBelosGMRES: has not been initialized.");
    Exceptions::amanzi_throw(msg);
  }

  Teuchos::RCP<Teuchos::ParameterList> pl =
    Teuchos::rcp(new Teuchos::ParameterList());
  pl->set("Num Blocks", krylov_dim_);
  pl->set("Maximum Iterations", max_itrs_);
  pl->set("Maximum Restarts", 2 * krylov_dim_ * max_itrs_);
  pl->set("Convergence Tolerance", tol_);
  pl->set("Output Frequency", 1);

  if (vo_->getVerbLevel() >= Teuchos::VERB_EXTREME) {
    pl->set("Verbosity",
            Belos::Errors + Belos::Warnings + Belos::IterationDetails +
              Belos::StatusTestDetails + Belos::Debug + Belos::TimingDetails);
  } else if (vo_->getVerbLevel() >= Teuchos::VERB_HIGH) {
    pl->set("Verbosity",
            Belos::Errors + Belos::Warnings + Belos::IterationDetails +
              Belos::StatusTestDetails);
  } else if (vo_->getVerbLevel() >= Teuchos::VERB_MEDIUM) {
    pl->set("Verbosity", Belos::Errors + Belos::Warnings);
  } else if (vo_->getVerbLevel() >= Teuchos::VERB_LOW) {
    pl->set("Verbosity", Belos::Errors + Belos::Warnings);
  }

  if (criteria_ & LIN_SOLVER_RELATIVE_RHS) {
    pl->set("Implicit Residual Scaling", "Norm of RHS");
    pl->set("Explicit Residual Scaling", "Norm of RHS");
  } else if (criteria_ & LIN_SOLVER_RELATIVE_RESIDUAL) {
    pl->set("Implicit Residual Scaling", "Norm of Initial Residual");
    pl->set("Explicit Residual Scaling", "Norm of Initial Residual");
  } else // if (criteria_ & LIN_SOLVER_ABSOLUTE_RESIDUAL)
  {
    pl->set("Implicit Residual Scaling", "None");
    pl->set("Explicit Residual Scaling", "None");
  }

  Teuchos::RCP<Belos::Operator<double>> Op =
    Teuchos::rcp(new AmanziBelosOp<Matrix, Vector>(this->m_));
  Teuchos::RCP<Belos::Operator<double>> Prec =
    Teuchos::rcp(new AmanziBelosOp<Matrix, Vector>(this->h_, true));
  Teuchos::RCP<Belos::MultiVec<double>> LHS =
    Teuchos::rcp(new CompositeMultiVector<Vector>(Teuchos::rcp(&hv, false)));
  Teuchos::RCP<const Belos::MultiVec<double>> RHS =
    Teuchos::rcp(new CompositeMultiVector<Vector>(
      Teuchos::rcp(const_cast<Vector*>(&v), false)));
  Teuchos::RCP<LinearProblem> problem =
    Teuchos::rcp(new LinearProblem(Op, LHS, RHS));

  problem->setLeftPrec(Prec);
  problem->setProblem();

  GmresSolver solver(problem, pl);
  Belos::ReturnType success = solver.solve();
  num_itrs_ = solver.getNumIters();
  residual_ = solver.achievedTol();

  if (success == Belos::Converged)
    returned_code_ = LIN_SOLVER_BELOS_SAYS_SUCCESS;
  else
    returned_code_ = LIN_SOLVER_BELOS_SAYS_FAIL;

  return returned_code_;
}

} // namespace AmanziSolvers
} // namespace Amanzi

#endif
