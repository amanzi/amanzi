/*
  Navier Stokes PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_NAVIER_STOKES_PK_HH_
#define AMANZI_NAVIER_STOKES_PK_HH_

#include <string>
#include <vector>

// TPLs
#include "Epetra_Vector.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

// Amanzi
#include "BDF1_TI.hh"
#include "EvaluatorPrimary.hh"
#include "Key.hh"
#include "PDE_Accumulation.hh"
#include "PDE_Abstract.hh"
#include "PDE_Elasticity.hh"
#include "PK.hh"
#include "PK_Factory.hh"
#include "PK_PhysicalBDF.hh"
#include "State.hh"
#include "TreeOperator.hh"
#include "TreeVector.hh"
#include "Units.hh"
#include "VerboseObject.hh"

#include "NavierStokesBoundaryFunction.hh"

namespace Amanzi {
namespace NavierStokes {

class NavierStokes_PK : public PK_PhysicalBDF {
 public:
  NavierStokes_PK(Teuchos::ParameterList& pk_tree,
                  const Teuchos::RCP<Teuchos::ParameterList>& glist,
                  const Teuchos::RCP<State>& S,
                  const Teuchos::RCP<TreeVector>& soln);

  NavierStokes_PK(const Teuchos::RCP<Teuchos::ParameterList>& glist,
                  const std::string& pk_list_name,
                  Teuchos::RCP<State> S,
                  const Teuchos::RCP<TreeVector>& soln);

  ~NavierStokes_PK() {};

  // methods required for PK interface
  virtual void Setup() final;
  virtual void Initialize() final;

  virtual double get_dt() final { return dt_; }
  virtual void set_dt(double dt) final { dt_ = dt; dt_desirable_ = dt_; }

  virtual bool AdvanceStep(double t_old, double t_new, bool reinit=false) final;
  virtual void CommitStep(double t_old, double t_new, const Tag& tag) final;
  virtual void CalculateDiagnostics(const Tag& tag) final {};

  virtual std::string name() { return passwd_; }

  // methods required for time integration interface
  // -- computes the non-linear functional f = f(t,u,udot) and related norm.
  void FunctionalResidual(const double t_old, double t_new, 
                  Teuchos::RCP<TreeVector> u_old, Teuchos::RCP<TreeVector> u_new, 
                  Teuchos::RCP<TreeVector> f);
  double ErrorNorm(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<const TreeVector> du);
  
  // -- management of the preconditioner
  int ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> pu);
  void UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> u, double dt);

  // -- check the admissibility of a solution
  //    override with the actual admissibility check
  bool IsAdmissible(Teuchos::RCP<const TreeVector> up) { return true; }

  // -- possibly modifies the predictor that is going to be used as a
  //    starting value for the nonlinear solve in the time integrator,
  //    the time integrator will pass the predictor that is computed
  //    using extrapolation and the time step that is used to compute
  //    this predictor this function returns true if the predictor was
  //    modified, false if not
  bool ModifyPredictor(double dt, 
                       Teuchos::RCP<const TreeVector> u0,
                       Teuchos::RCP<TreeVector> u) { return false; }

  // -- possibly modifies the correction, after the nonlinear solver (NKA)
  //    has computed it, will return true if it did change the correction,
  //    so that the nonlinear iteration can store the modified correction
  //    and pass it to NKA so that the NKA space can be updated
  AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
      ModifyCorrection(double dt,
                       Teuchos::RCP<const TreeVector> res,
                       Teuchos::RCP<const TreeVector> u,
                       Teuchos::RCP<TreeVector> du) {
    return AmanziSolvers::FnBaseDefs::CORRECTION_NOT_MODIFIED;
  }


  // -- calling this indicates that the time integration
  //    scheme is changing the value of the solution in state.
  void ChangedSolution() {
    pressure_eval_->SetChanged();
    fluid_velocity_eval_->SetChanged();
  }

  // other methods
  // --- management of boundary and source terms
  void ComputeOperatorBCs();

  // -- access
  Teuchos::RCP<BDF1_TI<TreeVector, TreeVectorSpace> > bdf1_dae() { return bdf1_dae_; }

 private:
  void UpdateSourceBoundaryData_(double t_old, double t_new);

 public:
  const Teuchos::RCP<Teuchos::ParameterList> glist_;
  Teuchos::RCP<Teuchos::ParameterList> ns_list_;
  Teuchos::RCP<const Teuchos::ParameterList> preconditioner_list_;
  Teuchos::RCP<const Teuchos::ParameterList> linear_solver_list_;
  Teuchos::RCP<Teuchos::ParameterList> ti_list_;

  int ncells_owned, ncells_wghost;
  int nfaces_owned, nfaces_wghost;
  int nnodes_owned, nnodes_wghost;

  double dt_, dt_next_, dt_desirable_;

 protected:
  // pointers to primary fields and their evaluators
  Teuchos::RCP<TreeVector> soln_;
  Teuchos::RCP<CompositeVector> soln_p_, soln_u_;

  Teuchos::RCP<EvaluatorPrimary<CompositeVector, CompositeVectorSpace> > pressure_eval_, fluid_velocity_eval_;

  // solvers
  Teuchos::RCP<Operators::TreeOperator> op_matrix_, op_preconditioner_, op_pc_solver_;
  Teuchos::RCP<Operators::PDE_Elasticity> op_matrix_elas_, op_preconditioner_elas_;
  Teuchos::RCP<Operators::PDE_Accumulation> op_matrix_acc_, op_preconditioner_acc_, op_mass_;
  Teuchos::RCP<Operators::PDE_Abstract> op_matrix_div_, op_matrix_grad_;
  Teuchos::RCP<Operators::PDE_Abstract> op_matrix_conv_, op_preconditioner_conv_;
  std::string solver_name_;
 
 private:
  std::string passwd_;
  int dim;

  Key pressure_key_, velocity_key_;

  // time integrators
  Teuchos::RCP<BDF1_TI<TreeVector, TreeVectorSpace> > bdf1_dae_;
  int num_itrs_;

  // boundary conditions
  std::vector<Teuchos::RCP<NavierStokesBoundaryFunction> > bcs_; 
  std::vector<Teuchos::RCP<Operators::BCs> > op_bcs_;

  // io
  Utils::Units units_;

  // factory registration
  static RegisteredPKFactory<NavierStokes_PK> reg_;
};

}  // namespace NavierStokes
}  // namespace Amanzi

#endif
