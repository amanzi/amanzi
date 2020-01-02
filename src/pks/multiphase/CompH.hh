/*
This is the multiphase flow component of the Amanzi code. 

Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided in the top-level COPYRIGHT file.

Authors: Quan Bui (mquanbui@math.umd.edu)
*/

#ifndef AMANZI_COMPH_PK_HH_
#define AMANZI_COMPH_PK_HH_

// TPLs
#include "Teuchos_RCP.hpp"

// Amanzi
#include "FlowBoundaryFunction.hh"
#include "PK_PhysicalBDF.hh"
#include "Tensor.hh"
#include "TI_Specs.hh"
#include "TreeVector.hh"
#include "State.hh"
#include "PDE_AdvectionUpwind.hh"
#include "PDE_Diffusion.hh"
#include "PDE_DiffusionFactory.hh"
#include "PDE_DiffusionFV.hh"
#include "PDE_DiffusionFVwithGravity.hh"
#include "PDE_Accumulation.hh"
#include "PK_DomainFunction.hh"
#include "PK_Factory.hh"
#include "UpwindFlux.hh"

// Amanzi::Multiphase
#include "CapillaryPressure.hh"
#include "MPCoeff.hh"
#include "MultiphaseTypeDefs.hh"
#include "WRMmp.hh"

namespace Amanzi {
namespace Multiphase {

class CompH_PK: public PK_PhysicalBDF {
 public:
  CompH_PK(Teuchos::ParameterList& pk_tree,
           const Teuchos::RCP<Teuchos::ParameterList>& global_list,
           const Teuchos::RCP<State>& S,
           const Teuchos::RCP<TreeVector>& soln);

  ~CompH_PK();

  // New interface for a PK
  virtual void Setup(const Teuchos::Ptr<State>& S) override {};
  virtual void Initialize(const Teuchos::Ptr<State>& S) override {
    InitializeFields();
    InitializeComponent();
    InitNextTI();
  }

  virtual double get_dt() override { return 1; }
  virtual void set_dt(double) override {};
  virtual bool AdvanceStep(double t_old, double t_new, bool reinit) override { return true; }
  virtual void CommitStep(double t_old, double t_new, const Teuchos::RCP<State>& S) override;
  virtual void CalculateDiagnostics(const Teuchos::RCP<State>& S) override {};
  virtual std::string name() override { return "hydrogen component"; }

  // Main methods of this PK
  void InitializeFields();
  void InitializeComponent();
  void InitNextTI();
  // void CommitState(const Teuchos::Ptr<State>& S);

  // Time integration interface new_mpc, implemented in Pressure_PK_TI.cc
  // computes the non-linear functional f = f(t,u,udot)
  virtual void FunctionalResidual(double t_old, double t_new, 
                                  Teuchos::RCP<TreeVector> u_old,
                                  Teuchos::RCP<TreeVector> u_new,
                                  Teuchos::RCP<TreeVector> f) override;

  // applies preconditioner to u and returns the result in Pu
  virtual int ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> pu) override;

  // updates the preconditioner
  virtual void UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h) override;

  // computes a norm on u-du and returns the result
  virtual double ErrorNorm(Teuchos::RCP<const TreeVector> u,
                           Teuchos::RCP<const TreeVector> du) override { return 0.0; }

  // check the admissibility of a solution
  // override with the actual admissibility check
  virtual bool IsAdmissible(Teuchos::RCP<const TreeVector> up) override { return true; }

  // possibly modifies the predictor that is going to be used as a
  // starting value for the nonlinear solve in the time integrator,
  // the time integrator will pass the predictor that is computed
  // using extrapolation and the time step that is used to compute
  // this predictor this function returns true if the predictor was
  // modified, false if not
  virtual bool ModifyPredictor(double h, Teuchos::RCP<const TreeVector> u0,
                               Teuchos::RCP<TreeVector> u) override { return false; }

  // possibly modifies the correction, after the nonlinear solver (NKA)
  // has computed it, will return true if it did change the correction,
  // so that the nonlinear iteration can store the modified correction
  // and pass it to NKA so that the NKA space can be updated
  virtual AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
  ModifyCorrection(double h, Teuchos::RCP<const TreeVector> res,
         Teuchos::RCP<const TreeVector> u,
         Teuchos::RCP<TreeVector> du) override {
    return AmanziSolvers::FnBaseDefs::CORRECTION_NOT_MODIFIED;
  }

  // experimental approach -- calling this indicates that the time
  // integration scheme is changing the value of the solution in
  // state.
  virtual void ChangedSolution() override {};

  void NumericalJacobian(double t_old, double t_new, Teuchos::RCP<const TreeVector> u, double eps);

  // methods to compute boundary and source terms
  //void SetAbsolutePermeabilityTensor(Teuchos::RCP<CompositeVector> Sw);
  void SetAbsolutePermeabilityTensor();
  void SetDiffusionTensor();
  void AddSourceTerms(CompositeVector& rhs) {};
  void ComputeBCs(bool stop);
  void ComputeBC_Pn(); 

  // access member functions
  Teuchos::RCP<Operators::PDE_AdvectionUpwind> OpPrec1() { return op_prec_pres_; }
  Teuchos::RCP<Operators::PDE_AdvectionUpwind> OpPrec2() { return op_prec_sat_; }
  Teuchos::RCP<Operators::PDE_AdvectionUpwind> OpPrec3() { return op_prec_rho_; }
  std::vector<Teuchos::RCP<Operators::PDE_Diffusion> >& Ops() { return ops_; }
  void SetJacobianType(std::string type) { jacobian_type_ = type; }

public:
  int ncells_owned_, ncells_wghost_;
  int nfaces_owned_, nfaces_wghost_;
  int ti_phase_counter_;
  int missed_bc_faces_, dirichlet_bc_faces_;

  Teuchos::ParameterList linear_operator_list_;
  Teuchos::ParameterList preconditioner_list_;
  Teuchos::RCP<Teuchos::ParameterList> comp_list_;
  Teuchos::ParameterList ti_list_;
  Teuchos::RCP<Teuchos::ParameterList> op_list_;
  Teuchos::ParameterList wrm_list_;

  double dT_, T_physics_, dTnext_;

private:
  // mesh structure and geometry
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  int dim_;
  double Cg_;

  // Stationary physical quantatities
  std::vector<WhetStone::Tensor> K_;
  std::vector<WhetStone::Tensor> D1_;
  AmanziGeometry::Point gravity_;
  double g_, atm_pressure_;
  double mu1_, mu2_, phi_, rho_;
  Teuchos::RCP<CompositeVector> rho_w_;

  // source and sink terms
  PK_DomainFunction* src_sink_;
  int src_sink_distribution_;

  // Time integration control
  TI_Specs* ti_specs_;
  TI_Specs ti_specs_generic_;

  std::string jacobian_type_, passwd_;

  // Verbose control
  VerboseObject* vo_;

  // boundary conditons
  std::vector<Teuchos::RCP<Flow::FlowBoundaryFunction> > bcs_; 
  Teuchos::RCP<Operators::BCs> op_bc_p_;
  Teuchos::RCP<Operators::BCs> op_bc_p_n_;
  Teuchos::RCP<Operators::BCs> op_bc_s_;
  Teuchos::RCP<Operators::BCs> op_bc_rhl_;
  
  // State and operators
  int error_control_, update_upwind_;
  double dT_desirable_;

  Teuchos::RCP<State> S_;

  Teuchos::RCP<MPCoeff> coef_w_;
  Teuchos::RCP<MPCoeff> coef_n_;
  Teuchos::RCP<CapillaryPressure> capillary_pressure_;

  Teuchos::RCP<Operators::PDE_DiffusionFVwithGravity> op1_matrix_;
  Teuchos::RCP<Operators::PDE_DiffusionFVwithGravity> op2_matrix_;
  Teuchos::RCP<Operators::PDE_Diffusion> op3_matrix_;
  Teuchos::RCP<Operators::PDE_Diffusion> op1_preconditioner_;
  Teuchos::RCP<Operators::PDE_Diffusion> op2_preconditioner_;
  Teuchos::RCP<Operators::PDE_Diffusion> op3_preconditioner_;
  Teuchos::RCP<Operators::PDE_AdvectionUpwind> op_prec_pres_;
  Teuchos::RCP<Operators::PDE_AdvectionUpwind> op_prec_rho_;
  Teuchos::RCP<Operators::PDE_AdvectionUpwind> op_prec_sat_;
  Teuchos::RCP<Operators::PDE_Accumulation> op_acc_;  
  Teuchos::RCP<Operators::PDE_Accumulation> op1_acc_;  
  Teuchos::RCP<Operators::PDE_Accumulation> op2_acc_;  

  // upwind operator
  Teuchos::RCP<Operators::UpwindFlux<MPCoeff> > upwind_w_;
  Teuchos::RCP<Operators::UpwindFlux<MPCoeff> > upwind_n_;

  std::vector<Teuchos::RCP<Operators::PDE_Diffusion> > ops_;
  typedef std::vector<Teuchos::RCP<Operators::PDE_Diffusion> >::iterator op_iter;
  typedef std::vector<Teuchos::RCP<Operators::Op> >::iterator local_op_iter;

  // The solution obtained from solving for pressure
  Teuchos::RCP<CompositeVector> sol_;
  //Teuchos::RCP<CompositeVector> saturation_phase1_;

  // solution tree vector
  Teuchos::RCP<TreeVector> soln_;

  // Darcy flux evaluator
  //Teuchos::RCP<PrimaryVariableFieldEvaluator> phase_velocity_eval_;

  // upwind flux
  Teuchos::RCP<CompositeVector> upwind_vw_;
  Teuchos::RCP<CompositeVector> upwind_vn_;
  Teuchos::RCP<CompositeVector> tmp_flux1_;
  Teuchos::RCP<CompositeVector> tmp_flux2_;

  //static RegisteredPKFactory<CompH_PK> reg_;

};

}  // namespase Flow
}  // namespace Amanzi

#endif

