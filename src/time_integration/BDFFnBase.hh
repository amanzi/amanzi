/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:

*/

//!

#ifndef AMANZI_BDFFNBASE_HH_
#define AMANZI_BDFFNBASE_HH_

#include "FnBaseDefs.hh"

namespace Amanzi {

// This is the interface definition for the BDF class.
// The nonlinear functional, preconditioner, and error
// functions must be derived from this class to be
// usable with BDF1_TI.

template <class Vector>
class BDFFnBase {
 public:
  virtual ~BDFFnBase() = default;

  // computes the non-linear functional f = f(t,u,udot)
  virtual void
  FunctionalResidual(double t_old, double t_new, Teuchos::RCP<Vector> u_old,
                     Teuchos::RCP<Vector> u_new, Teuchos::RCP<Vector> f) = 0;

  // applies preconditioner to u and returns the result in Pu
  virtual int ApplyPreconditioner(Teuchos::RCP<const Vector> u,
                                  Teuchos::RCP<Vector> Pu) = 0;

  // computes a norm on u-du and returns the result
  virtual double
  ErrorNorm(Teuchos::RCP<const Vector> u, Teuchos::RCP<const Vector> du) = 0;

  // updates the preconditioner
  virtual void
  UpdatePreconditioner(double t, Teuchos::RCP<const Vector> up, double h) = 0;

  // check the admissibility of a solution
  // override with the actual admissibility check
  virtual bool IsAdmissible(Teuchos::RCP<const Vector> up) = 0;

  // possibly modifies the predictor that is going to be used as a
  // starting value for the nonlinear solve in the time integrator,
  // the time integrator will pass the predictor that is computed
  // using extrapolation and the time step that is used to compute
  // this predictor this function returns true if the predictor was
  // modified, false if not
  virtual bool ModifyPredictor(double h, Teuchos::RCP<const Vector> u0,
                               Teuchos::RCP<Vector> u) = 0;

  // possibly modifies the correction, after the nonlinear solver (NKA)
  // has computed it, will return true if it did change the correction,
  // so that the nonlinear iteration can store the modified correction
  // and pass it to NKA so that the NKA space can be updated
  virtual AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
  ModifyCorrection(double h, Teuchos::RCP<const Vector> res,
                   Teuchos::RCP<const Vector> u, Teuchos::RCP<Vector> du) = 0;

  // update the continuation parameter
  virtual void UpdateContinuationParameter(double lambda){};

  // calling this indicates that the time
  // integration scheme is changing the value of the solution in
  // state.
  virtual void ChangedSolution() = 0;

  // experimental routine -- returns the number of linear iterations.
  virtual int ReportStatistics() { return 0; }
};

} // namespace Amanzi

#endif
