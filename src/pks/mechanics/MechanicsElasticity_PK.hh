/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*!

The conceptual PDE model for elasticity is

.. math::
  \frac{\partial (\rho \boldsymbol{u})}{\partial t} 
  + \boldsymbol{\nabla} \cdot (C \boldsymbol{u} \nabla \boldsymbol{u})
  =
  \rho \boldsymbol{g}

where 
:math:`\rho` is the fluid density [kg/m^3],
:math:`\boldsymbol{\sigma}` is the deviatoric stress tensor,
and
:math:`\boldsymbol{g}` is the gravity vector [:math:`m/s^2`].

*/

#ifndef AMANZI_MECHANICS_ELASTICITY_PK_HH_
#define AMANZI_MECHANICS_ELASTICITY_PK_HH_

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

#include "HydrostaticStressEvaluator.hh"
#include "MechanicsBoundaryFunction.hh"
#include "VolumetricStrainEvaluator.hh"

namespace Amanzi {
namespace Mechanics {

class MechanicsElasticity_PK : public PK_PhysicalBDF {
 public:
  MechanicsElasticity_PK(Teuchos::ParameterList& pk_tree,
                         const Teuchos::RCP<Teuchos::ParameterList>& glist,
                         const Teuchos::RCP<State>& S,
                         const Teuchos::RCP<TreeVector>& soln);

  ~MechanicsElasticity_PK(){};

  // methods required for PK interface
  virtual void Setup() final;
  virtual void Initialize() final;

  virtual double get_dt() final { return dt_; }
  virtual void set_dt(double dt) final { dt_ = dt; }

  virtual bool AdvanceStep(double t_old, double t_new, bool reinit = false) final;
  virtual void CommitStep(double t_old, double t_new, const Tag& tag) final;
  virtual void CalculateDiagnostics(const Tag& tag) final{};

  virtual std::string name() override { return passwd_; }

  // methods required for time integration interface
  // -- computes the non-linear functional f = f(t,u,udot) and related norm.
  void FunctionalResidual(const double t_old,
                          double t_new,
                          Teuchos::RCP<TreeVector> u_old,
                          Teuchos::RCP<TreeVector> u_new,
                          Teuchos::RCP<TreeVector> f) override;
  double ErrorNorm(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<const TreeVector> du) override;

  // -- management of the preconditioner
  int ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> pu) override;
  void UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> u, double dt) override;

  // -- check the admissibility of a solution
  //    override with the actual admissibility check
  bool IsAdmissible(Teuchos::RCP<const TreeVector> up) override { return true; }

  // -- possibly modifies the predictor that is going to be used as a
  //    starting value for the nonlinear solve in the time integrator,
  //    the time integrator will pass the predictor that is computed
  //    using extrapolation and the time step that is used to compute
  //    this predictor this function returns true if the predictor was
  //    modified, false if not
  bool ModifyPredictor(double dt, Teuchos::RCP<const TreeVector> u0, Teuchos::RCP<TreeVector> u) override { return false; }

  // -- possibly modifies the correction, after the nonlinear solver (NKA)
  //    has computed it, will return true if it did change the correction,
  //    so that the nonlinear iteration can store the modified correction
  //    and pass it to NKA so that the NKA space can be updated
  AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
  ModifyCorrection(double dt,
                   Teuchos::RCP<const TreeVector> res,
                   Teuchos::RCP<const TreeVector> u,
                   Teuchos::RCP<TreeVector> du) override
  {
    return AmanziSolvers::FnBaseDefs::CORRECTION_NOT_MODIFIED;
  }


  // -- calling this indicates that the time integration
  //    scheme is changing the value of the solution in state.
  void ChangedSolution() override { eval_->SetChanged(); }

  // other methods
  // --- management of boundary and source terms
  void ComputeOperatorBCs();

  // -- access
  Teuchos::RCP<BDF1_TI<TreeVector, TreeVectorSpace>> bdf1_dae() { return bdf1_dae_; }
  virtual Teuchos::RCP<Operators::Operator> my_operator(const Operators::OperatorType& type) override;
  virtual Teuchos::RCP<Operators::PDE_HelperDiscretization> my_pde(const Operators::PDEType& type) override { return op_matrix_elas_; }

 private:
  void UpdateSourceBoundaryData_(double t_old, double t_new);
  void AddGravityTerm_(CompositeVector& rhs);
  void AddPressureGradient_(CompositeVector& rhs);

 public:
  const Teuchos::RCP<Teuchos::ParameterList> glist_;
  Teuchos::RCP<Teuchos::ParameterList> ec_list_;
  Teuchos::RCP<const Teuchos::ParameterList> preconditioner_list_;
  Teuchos::RCP<const Teuchos::ParameterList> linear_solver_list_;
  Teuchos::RCP<Teuchos::ParameterList> ti_list_;

  int ncells_owned_, ncells_wghost_;
  int nfaces_owned_, nfaces_wghost_;
  int nnodes_owned_, nnodes_wghost_;

  double dt_, dt_next_, dt_desirable_;

  bool use_gravity_, biot_model_;

 protected:
  // pointers to primary fields and their evaluators
  Teuchos::RCP<TreeVector> soln_;
  Teuchos::RCP<CompositeVector> solution_;
  Teuchos::RCP<EvaluatorPrimary<CompositeVector, CompositeVectorSpace>> eval_;
  Teuchos::RCP<HydrostaticStressEvaluator> eval_hydro_stress_;
  Teuchos::RCP<VolumetricStrainEvaluator> eval_vol_strain_;

  // solvers
  Teuchos::RCP<Operators::Operator> op_matrix_, op_preconditioner_;
  Teuchos::RCP<Matrix<CompositeVector, CompositeVectorSpace>> op_pc_solver_;
  Teuchos::RCP<Operators::PDE_Elasticity> op_matrix_elas_, op_preconditioner_elas_;
  std::string solver_name_;

 private:
  std::string passwd_;
  int dim_;

  Key displacement_key_, hydrostatic_stress_key_;
  Key vol_strain_key_, prev_vol_strain_key_;
  Key young_modulus_key_, poisson_ratio_key_;
  Key particle_density_key_;

  // time integrators
  Teuchos::RCP<BDF1_TI<TreeVector, TreeVectorSpace>> bdf1_dae_;
  int num_itrs_;

  // boundary conditions
  std::vector<Teuchos::RCP<MechanicsBoundaryFunction>> bcs_;
  std::vector<Teuchos::RCP<Operators::BCs>> op_bcs_;
  int dirichlet_bc_;

  // io
  Utils::Units units_;

  // factory registration
  static RegisteredPKFactory<MechanicsElasticity_PK> reg_;
};

} // namespace Mechanics
} // namespace Amanzi

#endif
