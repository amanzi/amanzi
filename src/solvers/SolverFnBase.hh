/*
  This is the Nonlinear Solver component of the Amanzi code.

  Interface for an evaluator required for a nonlinear solver
  F(u) = 0.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/


#ifndef AMANZI_SOLVER_FN_BASE_
#define AMANZI_SOLVER_FN_BASE_

#include "Teuchos_RCP.hpp"

#include "FnBase_defs.hh"

namespace Amanzi {
namespace AmanziSolvers {

using namespace FnBaseDefs;


template<class Vector>
class SolverFnBase {
 public:
  // computes the non-linear functional r = F(u)
  virtual void Residual(const Teuchos::RCP<Vector>& u,
                        const Teuchos::RCP<Vector>& r) = 0;

  // preconditioner toolkit
  virtual void ApplyPreconditioner(const Teuchos::RCP<const Vector>& r,
                                   const Teuchos::RCP<Vector>& Pr) = 0;
  virtual void UpdatePreconditioner(const Teuchos::RCP<const Vector>& u) = 0;

  // error analysis
  virtual double ErrorNorm(const Teuchos::RCP<const Vector>& u,
                           const Teuchos::RCP<const Vector>& du) = 0;

  // allow PK to modify a correction
  virtual ModifyCorrectionResult ModifyCorrection(const Teuchos::RCP<const Vector>& r,
          const Teuchos::RCP<const Vector>& u,
          const Teuchos::RCP<Vector>& du) {
    return CORRECTION_NOT_MODIFIED;
  }

  // Check the admissibility of an inner iterate (ensures preconditions for
  // F(u) to be defined).
  virtual bool IsAdmissible(const Teuchos::RCP<const Vector>& up) {
    return true;
  }

  // bookkeeping for state
  virtual void ChangedSolution() {};

};

}  // namespace AmanziSolvers
}  // namespace Amanzi

#endif
