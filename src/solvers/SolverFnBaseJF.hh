/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  Solvers

  This decorator class wraps a nonlinear SolverFnBase with a class
  that replaces ApplyPreconditioner() with a Jacobian-free
  implementation of the inverse.

  Note this is a pass-through to the SolverFnBase in all but
  ApplyPreconditioner() and UpdatePreconditioner().
*/

#ifndef AMANZI_JF_SOLVER_FN_BASE_HH_
#define AMANZI_JF_SOLVER_FN_BASE_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "MatrixJF.hh"
#include "Matrix.hh"
#include "InverseFactory.hh"

#include "FnBaseDefs.hh"
#include "SolverFnBase.hh"

namespace Amanzi {
namespace AmanziSolvers {

template <class Vector, class VectorSpace>
class SolverFnBaseJF : public SolverFnBase<Vector> {
 public:
  SolverFnBaseJF(Teuchos::ParameterList& plist,
                 const Teuchos::RCP<SolverFnBase<Vector>> fn,
                 const VectorSpace& map);

  // -- Standard SolverFnBase interface.
  // computes the non-linear functional r = F(u)
  virtual void Residual(const Teuchos::RCP<Vector>& u, const Teuchos::RCP<Vector>& r)
  {
    fn_->Residual(u, r);
  }

  // preconditioner application
  virtual int
  ApplyPreconditioner(const Teuchos::RCP<const Vector>& r, const Teuchos::RCP<Vector>& Pr);

  // Update the preconditioner
  virtual void UpdatePreconditioner(const Teuchos::RCP<const Vector>& u0);

  // error norm
  virtual double
  ErrorNorm(const Teuchos::RCP<const Vector>& u, const Teuchos::RCP<const Vector>& du)
  {
    return fn_->ErrorNorm(u, du);
  }

  // Check the admissibility of an inner iterate (ensures preconditions for
  // F(u) to be defined).
  virtual bool IsAdmissible(const Teuchos::RCP<const Vector>& up) { return fn_->IsAdmissible(up); }

  // Hack a correction for some reason.
  virtual FnBaseDefs::ModifyCorrectionResult ModifyCorrection(const Teuchos::RCP<const Vector>& res,
                                                              const Teuchos::RCP<const Vector>& u,
                                                              const Teuchos::RCP<Vector>& du)
  {
    return fn_->ModifyCorrection(res, u, du);
  }

  // bookkeeping for state
  virtual void ChangedSolution() { fn_->ChangedSolution(); }

 protected:
  Teuchos::ParameterList plist_;
  Teuchos::RCP<MatrixJF<Vector, VectorSpace>> jf_mat_;
  Teuchos::RCP<Matrix<Vector, VectorSpace>> lin_op_;
  Teuchos::RCP<SolverFnBase<Vector>> fn_;

  double typical_u_;
};


// constructor
template <class Vector, class VectorSpace>
SolverFnBaseJF<Vector, VectorSpace>::SolverFnBaseJF(Teuchos::ParameterList& plist,
                                                    const Teuchos::RCP<SolverFnBase<Vector>> fn,
                                                    const VectorSpace& map)
  : plist_(plist), fn_(fn)
{
  typical_u_ = plist.get<double>("typical solution value", 1.0);

  // create the JF matrix, a Matrix<Vector,VectorSpace> class that
  // wraps the SolverFnBase as a Matrix doing residual evaluations to
  // approximate the action of the Jacobian.
  Teuchos::ParameterList jf_plist = plist_.sublist("JF matrix parameters");
  jf_mat_ = Teuchos::rcp(new MatrixJF<Vector, VectorSpace>(jf_plist, fn_, map));

  // Create the linear solver for that linear operator, which must be an
  // iterative method.
  Teuchos::ParameterList lin_plist = plist_.sublist("linear operator");
  if (!lin_plist.isParameter("iterative method")) {
    Errors::Message msg("JFNK \"linear operator\" sublist requires parameter \"iterative method\"");
    Exceptions::amanzi_throw(msg);
  }
  lin_op_ = AmanziSolvers::createIterativeMethod(lin_plist, jf_mat_);
}


// preconditioner application
template <class Vector, class VectorSpace>
int
SolverFnBaseJF<Vector, VectorSpace>::ApplyPreconditioner(const Teuchos::RCP<const Vector>& r,
                                                         const Teuchos::RCP<Vector>& Pr)
{
  return lin_op_->ApplyInverse(*r, *Pr);
}


// preconditioner update
template <class Vector, class VectorSpace>
void
SolverFnBaseJF<Vector, VectorSpace>::UpdatePreconditioner(const Teuchos::RCP<const Vector>& u0)
{
  fn_->UpdatePreconditioner(u0);
  jf_mat_->set_linearization_point(u0);
}

} // namespace AmanziSolvers
} // namespace Amanzi

#endif
