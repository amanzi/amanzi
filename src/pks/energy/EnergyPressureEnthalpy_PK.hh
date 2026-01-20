/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Energy PK

  Process kernel for energy equation in pressure-enthalpy variables.
*/

#ifndef AMANZI_ENERGY_PRESSURE_ENTHALPY_PK_HH_
#define AMANZI_ENERGY_PRESSURE_ENTHALPY_PK_HH_

// Amanzi
#include "BDF1_TI.hh"
#include "Energy_PK.hh"
#include "IEM.hh"
#include "ModelAssumptions.hh"
#include "PK_Factory.hh"
#include "StateArchive.hh"

// Energy

namespace Amanzi {
namespace Energy {

class EnergyPressureEnthalpy_PK : public Energy_PK {
 public:
  EnergyPressureEnthalpy_PK(Teuchos::ParameterList& pk_tree,
                            const Teuchos::RCP<Teuchos::ParameterList>& glist,
                            const Teuchos::RCP<State>& S,
                            const Teuchos::RCP<TreeVector>& soln);
  virtual ~EnergyPressureEnthalpy_PK() {};

  // methods required for PK interface
  virtual void Setup() final;
  virtual void Initialize() final;

  virtual bool AdvanceStep(double t_old, double t_new, bool reinit = false) final;
  virtual void CommitStep(double t_old, double t_new, const Tag& tag) final;
  virtual void FailStep(double t_old, double t_new, const Tag& tag) final;
  virtual void CalculateDiagnostics(const Tag& tag) final {};

  double get_dt() final { return dt_; }
  void set_dt(double dt) final { dt_ = dt; }

  virtual void FunctionalResidual(const double t_old,
                                  double t_new,
                                  Teuchos::RCP<const TreeVector> u_old,
                                  Teuchos::RCP<TreeVector> u_new,
                                  Teuchos::RCP<TreeVector> g) final;
  virtual double 
  ErrorNorm(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<const TreeVector> du) final;

  virtual void 
  UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double dt) override;

  // calling this indicates that the time integration
  // scheme is changing the value of the solution in state.
  void ChangedSolution() override { enthalpy_eval_->SetChanged(); }

  // additional method for the base class
  virtual void ComputeSecondaryBCs() override;

 private:
  void InitializeEnthalpy_();

 protected:
  const Teuchos::RCP<Teuchos::ParameterList> glist_;

  Key state_key_;
  Key bcs_flow_key_;
  std::string passwd_;
  Teuchos::RCP<EvaluatorPrimary<CompositeVector, CompositeVectorSpace>> enthalpy_eval_;

 private:
  // primary field
  const Teuchos::RCP<TreeVector> soln_;
  Teuchos::RCP<CompositeVector> solution;

  // operators and solvers
  Teuchos::RCP<Operators::PDE_Diffusion> op_matrix_diff_pres_;
  Teuchos::RCP<Operators::PDE_AdvectionUpwind> op_preconditioner_adv_enth_;

  // timestepping
  Teuchos::RCP<BDF1_TI<TreeVector, TreeVectorSpace>> bdf1_dae_;
  int num_itrs_;
  double dt_, dt_next_;

  Teuchos::RCP<StateArchive> archive_;

  static RegisteredPKFactory<EnergyPressureEnthalpy_PK> reg_;
};

} // namespace Energy
} // namespace Amanzi

#endif
