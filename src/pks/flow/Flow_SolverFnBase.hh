/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include "SolverFnBase.hh"

#include "Flow_PK.hh"

namespace Amanzi {
namespace Flow {

template<class Vector>
class Flow_SolverFnBase : public AmanziSolvers::SolverFnBase<Vector> {
 public:
  Flow_SolverFnBase(Teuchos::ParameterList& plist, Flow_PK* fn)
    : fn_(fn) {};

  virtual void
  Residual(const Teuchos::RCP<Vector>& u, const Teuchos::RCP<Vector>& r) {
    fn_->FunctionalResidual(0.0, 1.0e+98, u, u, r);
  }

  virtual int 
  ApplyPreconditioner(const Teuchos::RCP<const Vector>& u, const Teuchos::RCP<Vector>& hu) {
    return fn_->ApplyPreconditioner(u, hu);
  }

  virtual void
  UpdatePreconditioner(const Teuchos::RCP<const Vector>& u) {
    return fn_->UpdatePreconditioner(0.0, u, 1.0e+98);
  }

  virtual double
  ErrorNorm(const Teuchos::RCP<const Vector>& u, const Teuchos::RCP<const Vector>& du) {
    return fn_->ErrorNorm(u, du);
  }

  virtual void ChangedSolution() {
    fn_->ChangedSolution();
  }

 private:
  Flow_PK* fn_;
};

} // namespace Flow
} // namespace Amanzi

