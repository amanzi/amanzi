/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon (coonet@ornl.gov)
*/

//! Newton nonlinear solver

/*!

  Uses a classic Newton iteration.

  NOTE: this is only truely Newton if the provided SolverFnBase provides the
  true Jacobian in the ApplyJacobian() method.  If not, it is either inexact
  Newton or Picard, depending upon your preferences.

 */


#ifndef AMANZI_NEWTON_SOLVER_
#define AMANZI_NEWTON_SOLVER_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "SolverDefs.hh"
#include "SolverDefault.hh"

namespace Amanzi {
namespace AmanziSolvers {

template <class Vector, class VectorSpace>
class SolverNewton : public SolverDefault<Vector, VectorSpace> {
 public:
  SolverNewton(Teuchos::ParameterList& plist)
      : SolverDefault<Vector,VectorSpace>(plist)
  {
    // update the verbose options
    this->vo_ = Teuchos::rcp(new VerboseObject(name(), plist));
  }

  SolverNewton(Teuchos::ParameterList& plist,
               const Teuchos::RCP<SolverFnBase<Vector>>& fn,
               const VectorSpace& map)
      : SolverNewton(plist)
  {
    Init(fn, map);
  }

  virtual std::string name() const override { return "Solver::Newton"; }

 protected:
  //
  // Newton's ModifyCorrection() does basically nothing other than check that
  // the correction supplied by the predictor is admissible.
  virtual std::pair<MonitorStatus,double> ModifyCorrection_(Teuchos::RCP<Vector>& r,
          const Teuchos::RCP<Vector>& u, Teuchos::RCP<Vector>& du) override;

  
};

template <class Vector, class VectorSpace>
std::pair<MonitorStatus,double>
SolverNewton<Vector,VectorSpace>::ModifyCorrection_(Teuchos::RCP<Vector>& r,
        const Teuchos::RCP<Vector>& u, Teuchos::RCP<Vector>& du)
{
  // Hack the correction
  if (this->modify_correction_) {
    bool hacked = this->fn_->ModifyCorrection(r, u, du);
  }

  if (this->fn_->IsAdmissible(du)) return std::make_pair(MonitorStatus::CONTINUE, (double) -1.0);
  else return std::make_pair(MonitorStatus::INADMISSIBLE_SOLUTION, (double) -1.0);
}

} // namespace AmanziSolvers
} // namespace Amanzi

#endif
