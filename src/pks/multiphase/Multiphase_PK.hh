/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Quan Bui (mquanbui@math.umd.edu)
           Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*!

The conceptual PDE model for the isothermal multiphase flow include
transport equations for components and nonlinear algebraic constraints for the phase presence. 
At the moment we consider two phases (liquid and gas), multiple components, and one 
constraint.
Each transport equation has the following form:

.. math::
  \frac{\partial \Theta}{\partial t}
  + \nabla \cdot \boldsymbol{\Psi} = Q,

where 
:math:`\Theta` is the storage and 
:math:`\boldsymbol{\Psi}` is the flux.
The storage term sums up component amount across two phases, :math:`\alpha=l` for liquid
phase and :math:`\alpha=g` for gas phase:

.. math::
  \Theta = \phi \sum_\alpha \eta_\alpha x_\alpha s_\alpha

where
:math:`\phi` is porosity [-],
:math:`\eta` is molar density [mol/m^3],
:math:`x` is molar fraction of component [-], and
:math:`s` is phase saturation [-].

The flux includes advective and diffusion terms:

.. math::
  \boldsymbol{\Psi} 
  = -\sum_\alpha \eta_\alpha \left(\boldsymbol{q}_\alpha + D_\alpha \nabla x \right)

where
:math:`\boldsymbol{q}` is Darcy phase velocity,
:math:`D` is molecular duffusion coefficient.

The nonlinear algebraic constraint may have different forms. One of the available forms is

.. math::
  min (s_g, 1 - x_l - x_g) = 0.

It implies that if gas compounent is present then we must have :math:`x_l + x_g = 1`.

The PK provides three choices of primary variables. 
The first one includes pressure liquid, mole gas fraction, and saturation liquid.
The second one includes pressure liquid, molar gas density, and saturation liquid.
The third one is used for verification purposes and is based on the model in Jaffre's paper. 
This model describes two-phase two-component system with water and hydrogen. 

*/

#ifndef AMANZI_MULTIPHASE_PK_HH_
#define AMANZI_MULTIPHASE_PK_HH_

#include <set>
#include <vector>

// Amanzi
#include "BCs.hh"
#include "BDF1_TI.hh"
#include "BDFFnBase.hh"
#include "Evaluator.hh"
#include "EvaluatorPrimary.hh"
#include "FlattenedTreeOperator.hh"
#include "Key.hh"
#include "PDE_Accumulation.hh"
#include "PDE_AdvectionUpwindFactory.hh"
#include "PDE_DiffusionFactory.hh"
#include "PDE_DiffusionFVwithGravity.hh"
#include "PDE_DiffusionFVonManifolds.hh"
#include "PK_Factory.hh"
#include "PK_PhysicalBDF.hh"
#include "State.hh"
#include "TreeVector.hh"
#include "Upwind.hh"

// Multiphase
#include "EquationStructure.hh"
#include "MultiphaseEvaluator.hh"
#include "MultiphaseBoundaryFunction.hh"
#include "MultiphaseTypeDefs.hh"

namespace Amanzi {
namespace Multiphase {

class Multiphase_PK : public PK_PhysicalBDF {
 public:
  Multiphase_PK(Teuchos::ParameterList& pk_tree,
                const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                const Teuchos::RCP<State>& S,
                const Teuchos::RCP<TreeVector>& soln);

  ~Multiphase_PK(){};

  // method required for abstract PK interface
  virtual void Setup() override;
  virtual void Initialize() override;

  virtual double get_dt() override { return dt_; }
  virtual void set_dt(double dt) override { dt_ = dt; }

  virtual bool AdvanceStep(double t_old, double t_new, bool reinit) override;
  virtual void CommitStep(double t_old, double t_new, const Tag& tag) override;
  virtual void CalculateDiagnostics(const Tag& tag) override{};

  virtual std::string name() override { return "multiphase"; }

  // methods required for time integration interface
  // -- computes the non-linear functional f = f(t,u,udot) and related norm.
  virtual void FunctionalResidual(double t_old,
                                  double t_new,
                                  Teuchos::RCP<TreeVector> u_old,
                                  Teuchos::RCP<TreeVector> u_new,
                                  Teuchos::RCP<TreeVector> f) override;

  double ErrorNorm(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<const TreeVector> du) override;

  // -- preconditioner management
  virtual int
  ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> pu) override;
  virtual void
  UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double dt) override;

  // -- check the admissibility of a solution
  //    override with the actual admissibility check
  virtual bool IsAdmissible(Teuchos::RCP<const TreeVector> u) override { return false; }

  // -- possibly modifies the predictor that is going to be used as a
  //    starting value for the nonlinear solve in the time integrator,
  //    the time integrator will pass the predictor that is computed
  //    using extrapolation and the time step that is used to compute
  //    this predictor this function returns true if the predictor was
  //    modified, false if not
  virtual bool
  ModifyPredictor(double dt, Teuchos::RCP<const TreeVector> u0, Teuchos::RCP<TreeVector> u) override
  {
    return false;
  }

  // possibly modifies the correction, after the nonlinear solver (NKA)
  // has computed it, will return true if it did change the correction,
  // so that the nonlinear iteration can store the modified correction
  // and pass it to NKA so that the NKA space can be updated
  virtual AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
  ModifyCorrection(double h,
                   Teuchos::RCP<const TreeVector> res,
                   Teuchos::RCP<const TreeVector> u,
                   Teuchos::RCP<TreeVector> du) override;

  // calling this indicates that the time integration scheme is changing
  // the value of the solution in state.
  virtual void ChangedSolution() override;

  // multiphase submodels
  void PopulateBCs(int icomp, bool flag);
  void CheckCompatibilityBCs(const Key& keyr, const Key& gname);
  void ModifyEvaluators(int neqn);

  Teuchos::RCP<TreeVector> soln() { return soln_; }
  Teuchos::RCP<Operators::TreeOperator> op_tree_pc() { return op_preconditioner_; }

 protected:
  void InitializeFieldFromField_(const std::string& field0,
                                 const std::string& field1,
                                 bool call_evaluator);
  Teuchos::ParameterList MyRequire_(const Key& key, const std::string& owner);

 private:
  int InitMPSystem_(const std::string& eqn_name, int enq_id, int enq_num);
  void InitializeFields_();

  void PopulateSecondaryBCs_();

  Teuchos::RCP<CompositeVector> CreateCVforUpwind_();

  // debug tool
  WhetStone::DenseMatrix
  FiniteDifferenceJacobian_(double t_old,
                            double t_new,
                            Teuchos::RCP<const TreeVector> u_old,
                            Teuchos::RCP<const TreeVector> u_new,
                            double eps);

 protected:
  int ncells_owned_, ncells_wghost_;
  int nfaces_owned_, nfaces_wghost_;
  int dim_;

  std::string passwd_;
  double dt_, dt_next_;

  // unknowns
  Teuchos::RCP<TreeVector> soln_;
  std::vector<std::string> soln_names_;
  std::set<std::string> secondary_names_;
  std::vector<std::string> component_names_;
  int num_primary_, num_phases_;

  Teuchos::RCP<EvaluatorPrimary<CompositeVector, CompositeVectorSpace>> saturation_liquid_eval_;
  Teuchos::RCP<EvaluatorPrimary<CompositeVector, CompositeVectorSpace>> x_gas_eval_;

  // variable evaluators
  Teuchos::RCP<Evaluator> eval_tws_, eval_tcs_;

  // complimentarity problem
  std::string ncp_;

  // keys
  Key pressure_liquid_key_, x_liquid_key_, x_vapor_key_, x_gas_key_;
  Key saturation_liquid_key_, saturation_gas_key_, temperature_key_;
  Key energy_key_, prev_energy_key_;
  Key porosity_key_, pressure_gas_key_, pressure_vapor_key_;
  Key tcc_liquid_key_, tcc_gas_key_;
  Key permeability_key_, relperm_liquid_key_, relperm_gas_key_;
  Key advection_liquid_key_, advection_gas_key_;
  Key viscosity_liquid_key_, viscosity_gas_key_;
  Key vol_flowrate_liquid_key_, vol_flowrate_gas_key_;
  Key mol_density_liquid_key_, mol_density_gas_key_;
  Key mass_density_liquid_key_, mass_density_gas_key_;
  Key tws_key_, tcs_key_, prev_tws_key_, prev_tcs_key_;
  Key ncp_f_key_, ncp_g_key_, ncp_fg_key_;
  Key ie_rock_key_, ie_liquid_key_, ie_gas_key_;
  Key conductivity_key_, particle_density_key_;
  Key enthalpy_liquid_key_, enthalpy_gas_key_;
  Key advection_enthalpy_liquid_key_, advection_enthalpy_gas_key_;

  // matrix and preconditioner
  Teuchos::RCP<Operators::FlattenedTreeOperator> op_preconditioner_;
  Teuchos::RCP<Matrix<TreeVector, TreeVectorSpace>> op_pc_solver_;
  bool op_pc_assembled_;

  Teuchos::RCP<Operators::PDE_DiffusionFactory> fac_diffK_, fac_diffD_;
  Teuchos::RCP<Operators::PDE_AdvectionUpwindFactory> fac_adv_;
  std::vector<Teuchos::RCP<Operators::PDE_Diffusion>> pde_diff_K_;
  Teuchos::RCP<Operators::PDE_Diffusion> pde_diff_D_;

  std::map<std::string, bool> system_;
  std::vector<EquationStructure> eqns_;
  std::vector<std::vector<int>> eqns_flattened_;
  std::vector<std::string> eval_flattened_;

  // boundary conditions
  std::vector<Teuchos::RCP<MultiphaseBoundaryFunction>> bcs_;
  std::map<std::string, Teuchos::RCP<Operators::BCs>> op_bcs_;

  // physical parameters
  double mu_l_, mu_g_, rho_l_, eta_l_, mol_mass_H2O_;
  std::vector<WhetStone::Tensor> K_;
  std::vector<double> mol_diff_l_, mol_diff_g_, mol_mass_, kH_;

  // water retention models
  Teuchos::RCP<WRMmpPartition> wrm_;

  // upwind
  Teuchos::RCP<Operators::Upwind> upwind_, upwind_new_;

  // time integration
  std::vector<std::string> flux_names_;

  // io
  Utils::Units units_;
  Teuchos::RCP<Teuchos::ParameterList> mp_list_;

 private:
  int missed_bc_faces_;
  double smooth_mu_; // smoothing parameter
  bool flow_on_manifold_;

  // solvers and preconditioners
  bool cpr_enhanced_;

  Teuchos::RCP<const Teuchos::ParameterList> glist_;
  Teuchos::RCP<const Teuchos::ParameterList> pc_list_;
  Teuchos::RCP<const Teuchos::ParameterList> linear_operator_list_;
  Teuchos::RCP<Teuchos::ParameterList> ti_list_;

  // time integration
  int num_itrs_;
  Teuchos::RCP<BDF1_TI<TreeVector, TreeVectorSpace>> bdf1_dae_;

  // miscaleneous
  AmanziGeometry::Point gravity_;
  double g_;

 private:
  static RegisteredPKFactory<Multiphase_PK> reg_;
};

} // namespace Multiphase
} // namespace Amanzi
#endif
