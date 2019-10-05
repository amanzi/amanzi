/*
This is the multiphase flow component of the Amanzi code. 

Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided in the top-level COPYRIGHT file.

Authors: Quan Bui (mquanbui@math.umd.edu)
*/

#ifndef AMANZI_SATURATION_PK_HH_
#define AMANZI_SATURATION_PK_HH_

// TPLs
#include "Teuchos_RCP.hpp"

// Amanzi
#include "BDF1_TI.hh"
#include "FlowBoundaryFunction.hh"
#include "LinearOperatorFactory.hh"
#include "Operator.hh"
#include "PDE_Accumulation.hh"
#include "PDE_AdvectionUpwind.hh"
#include "PDE_Diffusion.hh"
#include "PDE_DiffusionFV.hh"
#include "PDE_DiffusionFactory.hh"
#include "PK_DomainFunction.hh"
#include "PK_Factory.hh"
#include "PK_PhysicalBDF.hh"
#include "primary_variable_field_evaluator.hh"
#include "State.hh"
#include "Tensor.hh"
#include "TI_Specs.hh"
#include "TreeVector.hh"
#include "Upwind.hh"
#include "UpwindFactory.hh"

// Amanzi::Multiphase
#include "RelativePermeability.hh"
#include "WaterRetentionModel.hh"
#include "FractionalFlow.hh"
#include "CapillaryPressure.hh"

namespace Amanzi {
namespace Multiphase {

class Saturation_PK: public PK_PhysicalBDF {
public:
  Saturation_PK(Teuchos::ParameterList& pk_tree,
                const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                const Teuchos::RCP<State>& S,
                const Teuchos::RCP<TreeVector>& soln);

  ~Saturation_PK();

  // New interface for a PK
  virtual void Setup();
  virtual void Initialize() {
    InitializeFields();
    InitializeSaturation();
    InitTimeInterval(*ti_list_);
  }

  virtual double get_dt() {return 1.0;}
  virtual void set_dt(double) {};
  virtual bool AdvanceStep(double t_old, double t_new, bool reinit);
  virtual void CommitStep(double t_old, double t_new);
  virtual void CalculateDiagnostics() {};
  virtual std::string name(){return "saturation";}

  // Main methods of this PK
  void InitializeFields();
  void InitializeSaturation();
  void InitTimeInterval(Teuchos::ParameterList& ti_list);
  void InitNextTI(double T0, double dT0, TI_Specs& ti_specs);
  void CommitState(const Teuchos::Ptr<State>& S);

  bool Advance(double dT_MPC, double& dT_actual);

  // Time integration interface new_mpc, implemented in Saturation_PK_TI.cc
  // computes the non-linear functional f = f(t,u,udot)
  virtual void Functional(double t_old, double t_new, 
                          Teuchos::RCP<TreeVector> u_old,
                          Teuchos::RCP<TreeVector> u_new,
                          Teuchos::RCP<TreeVector> f);

  // applies preconditioner to u and returns the result in Pu
  virtual int ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, 
                                  Teuchos::RCP<TreeVector> Pu);

  // updates the preconditioner
  virtual void UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h);

  // computes a norm on u-du and returns the result
  virtual double ErrorNorm(Teuchos::RCP<const TreeVector> u,
                           Teuchos::RCP<const TreeVector> du);

  // check the admissibility of a solution
  // override with the actual admissibility check
  virtual bool IsAdmissible(Teuchos::RCP<const TreeVector> up) { return true; }

  // possibly modifies the predictor that is going to be used as a
  // starting value for the nonlinear solve in the time integrator,
  // the time integrator will pass the predictor that is computed
  // using extrapolation and the time step that is used to compute
  // this predictor this function returns true if the predictor was
  // modified, false if not
  virtual bool ModifyPredictor(double h, Teuchos::RCP<const TreeVector> u0,
                               Teuchos::RCP<TreeVector> u) { return false; }

  // possibly modifies the correction, after the nonlinear solver (NKA)
  // has computed it, will return true if it did change the correction,
  // so that the nonlinear iteration can store the modified correction
  // and pass it to NKA so that the NKA space can be updated
  virtual AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
  ModifyCorrection(double h, Teuchos::RCP<const TreeVector> res,
                   Teuchos::RCP<const TreeVector> u,
                   Teuchos::RCP<TreeVector> du);

  // experimental approach -- calling this indicates that the time
  // integration scheme is changing the value of the solution in
  // state.
  virtual void ChangedSolution() {};

  virtual int ReportStatistics() {return ls_itrs_;}

  void LineSearch(Teuchos::RCP<CompositeVector> ds);

  void ClipSaturation(Teuchos::RCP<CompositeVector> s, double tol);

  // misc stuff
  void ResetPKtimes(double T0, double dT0) { T_physics = T0; dT = dT0; }
  void AddSourceTerm(CompositeVector& rhs, const CompositeVector& fractional_flow);
  void ComputeBCs();

  // IO methods
  void ProcessParameterList(Teuchos::ParameterList& list);
  void ProcessSublistTimeIntegration(Teuchos::ParameterList& list, const std::string name, TI_Specs& ti_specs);
  virtual void AddSourceTerms(CompositeVector& rhs);
  void ProcessStringSourceDistribution(const std::string name, int* method);
  void ProcessStringTimeIntegration(const std::string name, int* method);
  void ProcessStringErrorOptions(Teuchos::ParameterList& list, int* control);
  void ProcessSublistTimeInterval(Teuchos::ParameterList& ti_list,  TI_Specs& ti_specs);
  std::string FindStringLinearSolver(const Teuchos::ParameterList& plist);
  std::string FindStringPreconditioner(const Teuchos::ParameterList& list);
  void OutputTimeHistory(const Teuchos::ParameterList& plist, std::vector<dt_tuple>& dT_history);

  void SetAbsolutePermeabilityTensor();
  void ComputeBC_Pc();

  void AnalyticJacobian(double Tp, Teuchos::RCP<const TreeVector> u, double dTp);
  void NumericalJacobian(double Tp, Teuchos::RCP<const TreeVector> u, double dTp);

  // access methods
  Teuchos::RCP<Operators::Operator> op_prec() { return op_preconditioner_->global_operator(); }

 public:
  int ncells_owned, ncells_wghost;
  int nfaces_owned, nfaces_wghost;

  double dT, T_physics, dTnext;
  std::vector<WhetStone::Tensor> K_;

  Teuchos::RCP<Teuchos::ParameterList> linear_operator_list_;
  Teuchos::RCP<Teuchos::ParameterList> preconditioner_list_;
  Teuchos::RCP<Teuchos::ParameterList> ti_list_;
  Teuchos::RCP<Teuchos::ParameterList> mp_list_, wrm_list_;

  int ti_phase_counter;
  int missed_bc_faces_, dirichlet_bc_faces_;

 private:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  int dim_;
  bool include_capillary_;

  // Verbose control
  VerboseObject* vo_;

  // boundary conditons
  std::vector<int> bc_submodel;
  Flow::FlowBoundaryFunction* bc_saturation;
  std::vector<Teuchos::RCP<Flow::FlowBoundaryFunction> > bcs_; 

  // time integration phases
  TI_Specs ti_specs_generic_;
  TI_Specs ti_specs_igs_;
  TI_Specs ti_specs_sss_;
  TI_Specs ti_specs_trs_;
  TI_Specs* ti_specs_;

  // source and sink terms
  PK_DomainFunction* src_sink_;
  int src_sink_distribution_;

  // Stationary physical quantatities
  AmanziGeometry::Point gravity_;
  std::vector<double> rho_, mu_;
  double g_, atm_pressure_;
  double phi_;

  std::string passwd_, jacobian_type_;

  int error_control_, update_upwind_, ls_itrs_;
  double dT_desirable_;

  // State, vectors, and operators
  Teuchos::RCP<State> S_;

  Teuchos::RCP<CompositeVector> tot_mobility_;
  Teuchos::RCP<CompositeVector> fractional_flow_;
  Teuchos::RCP<CompositeVector> dfw_dS_;
  Teuchos::RCP<CompositeVector> darcy_flux_;

  Teuchos::RCP<RelativePermeability> rel_perm_w_;
  Teuchos::RCP<RelativePermeability> rel_perm_n_;
  Teuchos::RCP<FractionalFlow> frac_flow_;
  Teuchos::RCP<CapillaryPressure> capillary_pressure_;
  //Teuchos::RCP<Operators::Operator> op_;
  Teuchos::RCP<Operators::PDE_Accumulation> op_acc_;
  Teuchos::RCP<Operators::PDE_AdvectionUpwind> op_matrix_;
  Teuchos::RCP<Operators::PDE_Diffusion> op1_matrix_;
  Teuchos::RCP<Operators::PDE_AdvectionUpwind> op_preconditioner_;
  Teuchos::RCP<Operators::PDE_AdvectionUpwind> op_sum_;
  Teuchos::RCP<Operators::PDE_AdvectionUpwind> op_sum1_;
  Teuchos::RCP<Operators::PDE_Diffusion> op1_preconditioner_;
  Teuchos::RCP<Operators::PDE_DiffusionFV> op_pres_pc_;
  Teuchos::RCP<Operators::BCs> op_bc_;
  Teuchos::RCP<Operators::BCs> op_bc_pc_;
  Teuchos::RCP<Operators::BCs> op_bc_tmp_;
  Teuchos::RCP<Operators::Upwind<FractionalFlow> > upwind_;
  Teuchos::RCP<Operators::Upwind<RelativePermeability> > upwind_n_;
  Teuchos::RCP<Operators::Upwind<CapillaryPressure> > upwind_pc_;

  Teuchos::RCP<BDF1_TI<TreeVector, TreeVectorSpace> > bdf1_dae;  // Time integrators

  // The solution obtained from solving for pressure
  Teuchos::RCP<CompositeVector> sol_;

  // solution tree vector, required for the mpc interface
  // I don't know what this will do yet
  Teuchos::RCP<TreeVector> soln_;

  AmanziSolvers::LinearOperatorFactory<Operators::Operator, CompositeVector, CompositeVectorSpace> sfactory;
  Teuchos::RCP<AmanziSolvers::LinearOperator<Operators::Operator, CompositeVector, CompositeVectorSpace> >
     solver;

  static RegisteredPKFactory<Saturation_PK> reg_;
};

}  // namespase Flow
}  // namespace Amanzi

#endif

