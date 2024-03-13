/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! A very simple nonlinear continuation method.
/*!

Continuation methods are useful when the nonlinearity can be controlled by a
single simple parameter.  In this method, the nonlinear problem is solved with
a less-nonlinear value of the parameter, and the solution of that is used as
the initial guess to solve a harder problem.  As each successive problem is
solved, the continuation parameter is changed closer and closer to the true
value.

Few if any PKs support this method currently -- it requires the PK to provide more
interface about how to update the continuation parameter.

.. _solver-continuation-spec:
.. admonition:: solver-continuation-spec

    * `"nonlinear tolerance`" ``[double]`` **1.e-6** defines the required error
      tolerance. The error is calculated by a PK.

    * `"number of continuation steps`" ``[int]`` **5** How many steps to take
      from initial parameter to final parameter.

    * `"inner solver`" ``[solver-typed-spec]`` A Solver_, used at each step.

*/


#ifndef AMANZI_CONTINUATION_SOLVER_
#define AMANZI_CONTINUATION_SOLVER_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "VerboseObject.hh"

#include "Solver.hh"
#include "SolverFactory.hh"
#include "SolverFnBase.hh"
#include "SolverDefs.hh"

namespace Amanzi {
namespace AmanziSolvers {

template <class Vector, class VectorSpace>
class SolverContinuation : public Solver<Vector, VectorSpace> {
 public:
  SolverContinuation(Teuchos::ParameterList& plist) : plist_(plist){};

  SolverContinuation(Teuchos::ParameterList& plist,
                     const Teuchos::RCP<SolverFnBase<Vector>>& fn,
                     const VectorSpace& map)
    : plist_(plist)
  {
    Init(fn, map);
  }

  void Init(const Teuchos::RCP<SolverFnBase<Vector>>& fn, const VectorSpace& map);

  int Solve(const Teuchos::RCP<Vector>& u)
  {
    returned_code_ = Solve_(u);
    return (returned_code_ >= 0) ? 0 : 1;
  }

  // mutators
  void set_tolerance(double tol)
  {
    tol_ = tol;
    solver_->set_tolerance(tol);
  }
  void set_pc_lag(int pc_lag) { solver_->set_pc_lag(pc_lag); }

  // access
  double tolerance() { return tol_; }
  double residual() { return solver_->residual(); }
  int num_itrs() { return num_itrs_; }
  int pc_calls() { return 0; }
  int pc_updates() { return 0; }
  int returned_code() { return returned_code_; }
  std::vector<std::pair<double, double>>& history() { return solver_->history(); }

 private:
  void Init_();
  int Solve_(const Teuchos::RCP<Vector>& u);

 protected:
  Teuchos::ParameterList plist_;
  Teuchos::RCP<SolverFnBase<Vector>> fn_;
  Teuchos::RCP<Solver<Vector, VectorSpace>> solver_;

  Teuchos::RCP<VerboseObject> vo_;

 private:
  double tol_;
  int num_itrs_;
  double n_cont_steps_;
  int returned_code_;
};


/* ******************************************************************
* Public Init method.
****************************************************************** */
template <class Vector, class VectorSpace>
void
SolverContinuation<Vector, VectorSpace>::Init(const Teuchos::RCP<SolverFnBase<Vector>>& fn,
                                              const VectorSpace& map)
{
  fn_ = fn;
  Init_();
  solver_->Init(fn, map);
}


/* ******************************************************************
* Initialization of the NKA solver.
****************************************************************** */
template <class Vector, class VectorSpace>
void
SolverContinuation<Vector, VectorSpace>::Init_()
{
  tol_ = plist_.get<double>("nonlinear tolerance", 1.e-6);
  n_cont_steps_ = plist_.get<int>("number of continuation steps", 5);

  SolverFactory<Vector, VectorSpace> fac;
  solver_ = fac.Create(plist_.sublist("inner solver"));

  // update the verbose options
  vo_ = Teuchos::rcp(new VerboseObject("Solver::Continuation", plist_));
}


/* ******************************************************************
* The body of NKA solver
****************************************************************** */
template <class Vector, class VectorSpace>
int
SolverContinuation<Vector, VectorSpace>::Solve_(const Teuchos::RCP<Vector>& u)
{
  Teuchos::OSTab tab = vo_->getOSTab();

  // initialize the iteration counter
  num_itrs_ = 0;

  int itr = 0;
  do {
    if (vo_->os_OK(Teuchos::VERB_MEDIUM)) {
      *vo_->os() << "taking continuation step " << itr << " of " << n_cont_steps_ << std::endl;
    }

    // continuation parameter ranges from 1 (smooth problem) to 0 (true problem)
    double lambda = (n_cont_steps_ - itr) / n_cont_steps_;
    fn_->UpdateContinuationParameter(lambda);

    // solve
    int ierr = solver_->Solve(u);

    num_itrs_ += solver_->num_itrs();

    if (ierr) { return solver_->returned_code(); }

    itr++;
  } while (itr <= n_cont_steps_);
  return itr;
}


} // namespace AmanziSolvers
} // namespace Amanzi

#endif
