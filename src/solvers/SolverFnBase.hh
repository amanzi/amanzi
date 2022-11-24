/*
  Solvers

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)

  Interface for an evaluator required for a nonlinear solver
  F(u) = 0.
*/


#ifndef AMANZI_SOLVER_FN_BASE_
#define AMANZI_SOLVER_FN_BASE_

#include "Teuchos_RCP.hpp"

#include "FnBaseDefs.hh"

namespace Amanzi {
namespace AmanziSolvers {

template <class Vector>
class SolverFnBase {
 public:
  virtual ~SolverFnBase() = default;

  // computes the non-linear functional r = F(u)
  virtual void Residual(const Teuchos::RCP<Vector>& u, const Teuchos::RCP<Vector>& r) = 0;

  // preconditioner toolkit
  virtual int
  ApplyPreconditioner(const Teuchos::RCP<const Vector>& r, const Teuchos::RCP<Vector>& Pr) = 0;
  virtual void UpdatePreconditioner(const Teuchos::RCP<const Vector>& u) = 0;

  // error analysis
  virtual double
  ErrorNorm(const Teuchos::RCP<const Vector>& u, const Teuchos::RCP<const Vector>& du) = 0;

  // allow PK to modify a correction
  virtual FnBaseDefs::ModifyCorrectionResult ModifyCorrection(const Teuchos::RCP<const Vector>& r,
                                                              const Teuchos::RCP<const Vector>& u,
                                                              const Teuchos::RCP<Vector>& du)
  {
    return FnBaseDefs::CORRECTION_NOT_MODIFIED;
  }

  virtual void UpdateContinuationParameter(double lambda){};

  // Check the admissibility of an inner iterate (ensures preconditions for
  // F(u) to be defined).
  virtual bool IsAdmissible(const Teuchos::RCP<const Vector>& up) { return true; }

  // bookkeeping for state
  virtual void ChangedSolution(){};
};

} // namespace AmanziSolvers
} // namespace Amanzi

#endif
