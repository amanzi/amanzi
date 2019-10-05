/*
This is the multiphase flow component of the Amanzi code. 

Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided in the top-level COPYRIGHT file.

Authors: Quan Bui (mquanbui@math.umd.edu)
*/

#ifndef AMANZI_COMP_PK_HH_
#define AMANZI_COMP_PK_HH_

// Trilinos include
#include "Teuchos_RCP.hpp"

// Basic data structure include
#include "Tensor.hh"
#include "TreeVector.hh"

// Time integration include
#include "FnTimeIntegratorPK.hh"
#include "TI_Specs.hh"

// General include
#include "State.hh"
#include "OperatorDiffusionFactory.hh"
#include "OperatorAdvection.hh"
#include "OperatorDiffusionFV.hh"
#include "PK_Factory.hh"
#include "UpwindFlux.hh"
#include "OperatorDiffusion.hh"
#include "OperatorAccumulation.hh"

// Specific include for this PK
#include "RelativePermeability.hh"
#include "CapillaryPressure.hh"
#include "MultiphaseTypeDefs.hh"
#include "WaterRetentionModel.hh"

#include "FlowBoundaryFunction.hh"
#include "FlowDomainFunction.hh"
#include "Flow_BC_Factory.hh"
#include "Flow_SourceFactory.hh"

namespace Amanzi {
namespace Multiphase {
//class State;

class Comp_PK: public FnTimeIntegratorPK {
public:
  Comp_PK(Teuchos::ParameterList& pk_tree,
                    const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                    const Teuchos::RCP<State>& S,
                    const Teuchos::RCP<TreeVector>& soln);

  ~Comp_PK();

  // New interface for a PK
  virtual void Setup(){};
  virtual void Initialize() {
      InitializeFields();
      InitializeComponent();
      InitNextTI();
  }

  virtual double get_dt(){return 1;}
  virtual void set_dt(double){};
  virtual bool AdvanceStep(double t_old, double t_new, bool reinit){};
  virtual void CommitStep(double t_old, double t_new);
  virtual void CalculateDiagnostics(){};
  virtual std::string name(){return comp_name_;}

  // Main methods of this PK
  void InitializeFields();
  void InitializeComponent();
  void InitNextTI();
  int Comp_ID() {return comp_id_;}
  //void CommitState(const Teuchos::Ptr<State>& S);

  // Time integration interface new_mpc, implemented in Pressure_PK_TI.cc
  // computes the non-linear functional f = f(t,u,udot)
  virtual void Functional(double t_old, double t_new, 
                          Teuchos::RCP<TreeVector> u_old,
                      Teuchos::RCP<TreeVector> u_new,
                          Teuchos::RCP<TreeVector> f);

  // applies preconditioner to u and returns the result in Pu
  virtual int ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, 
                                   Teuchos::RCP<TreeVector> Pu);

  // updates the preconditioner
  virtual void UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up,
                    double h);

  // computes a norm on u-du and returns the result
  virtual double ErrorNorm(Teuchos::RCP<const TreeVector> u,
           Teuchos::RCP<const TreeVector> du) {
  }

  // check the admissibility of a solution
  // override with the actual admissibility check
  virtual bool IsAdmissible(Teuchos::RCP<const TreeVector> up) {
  }

  // possibly modifies the predictor that is going to be used as a
  // starting value for the nonlinear solve in the time integrator,
  // the time integrator will pass the predictor that is computed
  // using extrapolation and the time step that is used to compute
  // this predictor this function returns true if the predictor was
  // modified, false if not
  virtual bool ModifyPredictor(double h, Teuchos::RCP<const TreeVector> u0,
         Teuchos::RCP<TreeVector> u) {
  }

  // possibly modifies the correction, after the nonlinear solver (NKA)
  // has computed it, will return true if it did change the correction,
  // so that the nonlinear iteration can store the modified correction
  // and pass it to NKA so that the NKA space can be updated
  virtual AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
  ModifyCorrection(double h, Teuchos::RCP<const TreeVector> res,
         Teuchos::RCP<const TreeVector> u,
         Teuchos::RCP<TreeVector> du) {
  }

  // experimental approach -- calling this indicates that the time
  // integration scheme is changing the value of the solution in
  // state.
  virtual void ChangedSolution() {}

  void NumericalJacobian(double t_old, double t_new, Teuchos::RCP<const TreeVector> u, double eps);

  // methods to compute boundary and source terms
  //void SetAbsolutePermeabilityTensor(Teuchos::RCP<CompositeVector> Sw);
  void SetAbsolutePermeabilityTensor();
  void SetDiffusionTensor();
  void AddSourceTerms(CompositeVector& rhs) {};
  void ComputeBCs();
  void ComputeBC_Pc(); 
  void DeriveFaceValuesFromCellValues(const Epetra_MultiVector& ucells, Epetra_MultiVector& ufaces,
                                      const std::vector<int>& bc_model, const std::vector<double>& bc_value);

  // io members, implemented in Pressure_PK_IO.cc
  void ProcessParameterList(Teuchos::ParameterList& list);
  void ProcessSublistTimeIntegration(Teuchos::ParameterList& list, const std::string name, TI_Specs& ti_specs);
  void ProcessStringSourceDistribution(const std::string name, int* method);
  void ProcessStringTimeIntegration(const std::string name, int* method);
  void ProcessStringErrorOptions(Teuchos::ParameterList& list, int* control);
  void ProcessSublistTimeInterval(Teuchos::ParameterList& ti_list,  TI_Specs& ti_specs);
  std::string FindStringLinearSolver(const Teuchos::ParameterList& plist);
  std::string FindStringPreconditioner(const Teuchos::ParameterList& list);
  void OutputTimeHistory(const Teuchos::ParameterList& plist, std::vector<dt_tuple>& dT_history);

  // access member functions
  Teuchos::RCP<Operators::PDE_DiffusionFV> op_pres_prec() { return op_pres_prec_; }
  Teuchos::RCP<Operators::PDE_DiffusionFV> op_fuga_prec() { return op4_preconditioner_; }
  Teuchos::RCP<Operators::PDE_AdvectionUpwind> op_sat_prec() { return op_sat_prec_; }
  std::vector<Teuchos::RCP<Operators::PDE_Diffusion> >& Ops() { return ops_; }

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
  Teuchos::ParameterList eos_list_;
  Teuchos::ParameterList wrm_list_;

  double dT_, T_physics_, dTnext_;

private:
  // mesh structure and geometry
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  int dim_, comp_id_;
  double comp_coeff_;

  // Stationary physical quantatities
  std::vector<WhetStone::Tensor> K_;
  std::vector<WhetStone::Tensor> D1_;
  std::vector<WhetStone::Tensor> D2_;
  //Teuchos::RCP<Epetra_Vector> Kxy_;
  AmanziGeometry::Point gravity_;
  double g_, atm_pressure_;
  double rho_, mu_, phi_;

  // source and sink terms
  Flow::FlowDomainFunction* src_sink_;
  int src_sink_distribution_;

  // Time integration control
  TI_Specs* ti_specs_;
  TI_Specs ti_specs_generic_;

  std::string passwd_, comp_name_;
  std::string molar_fraction_name1_, molar_fraction_name2_;
  std::string diff_tensor_name1_, diff_tensor_name2_;

  // Verbose control
  VerboseObject* vo_;

  // boundary conditons
  std::vector<int> bc_model_p_, bc_submodel_p_;
  std::vector<double> bc_value_p_, bc_coef_p_;
  std::vector<int> bc_model_pc_, bc_submodel_pc_;
  std::vector<double> bc_value_pc_, bc_coef_pc_;
  std::vector<int> bc_model_s_, bc_submodel_s_;
  std::vector<double> bc_value_s_, bc_coef_s_;
  std::vector<int> bc_model_cfw_, bc_submodel_cfw_;
  std::vector<double> bc_value_cfw_, bc_coef_cfw_;
  std::vector<int> bc_model_cfn_, bc_submodel_cfn_;
  std::vector<double> bc_value_cfn_, bc_coef_cfn_;
  std::vector<int> bc_model_xw_, bc_submodel_xw_;
  std::vector<double> bc_value_xw_, bc_coef_xw_;
  std::vector<int> bc_model_xn_, bc_submodel_xn_;
  std::vector<double> bc_value_xn_, bc_coef_xn_;
  std::vector<double> bc_mixed_;

  Flow::FlowBoundaryFunction* bc_pressure_;
  Flow::FlowBoundaryFunction* bc_saturation_;
  //Flow::FlowBoundaryFunction* bc_flux_phase1_;
  //Flow::FlowBoundaryFunction* bc_flux_phase2_;
  Flow::FlowBoundaryFunction* bc_mfw_;
  Flow::FlowBoundaryFunction* bc_mfn_;
  
  // State and operators
  int error_control_, update_upwind_;
  double dT_desirable_;

  Teuchos::RCP<State> S_;

  Teuchos::RCP<RelativePermeability> rel_perm_w_;
  Teuchos::RCP<RelativePermeability> rel_perm_n_;
  Teuchos::RCP<CapillaryPressure> capillary_pressure_;
  Teuchos::RCP<CapillaryPressure> capillary_pressure_old_;

  Teuchos::RCP<Operators::PDE_Diffusion> op1_matrix_;
  Teuchos::RCP<Operators::PDE_Diffusion> op2_matrix_;
  Teuchos::RCP<Operators::PDE_Diffusion> op3_matrix_;
  Teuchos::RCP<Operators::PDE_Diffusion> op4_matrix_;
  Teuchos::RCP<Operators::PDE_Diffusion> op1_preconditioner_;
  Teuchos::RCP<Operators::PDE_AdvectionUpwind> op2_preconditioner_;
  Teuchos::RCP<Operators::PDE_AdvectionUpwind> op5_preconditioner_;
  Teuchos::RCP<Operators::PDE_AdvectionUpwind> op6_preconditioner_;
  Teuchos::RCP<Operators::PDE_Diffusion> op3_preconditioner_;
  Teuchos::RCP<Operators::PDE_DiffusionFV> op4_preconditioner_;
  Teuchos::RCP<Operators::PDE_DiffusionFV> op_pres_prec_;
  Teuchos::RCP<Operators::PDE_AdvectionUpwind> op_sat_prec_;
  Teuchos::RCP<Operators::PDE_Accumulation> op_pres_acc_;
  Teuchos::RCP<Operators::PDE_Accumulation> op_sat_acc_;
  Teuchos::RCP<Operators::PDE_Accumulation> op_fug_acc_;
  Teuchos::RCP<Operators::BCs> op_bc_s_;
  Teuchos::RCP<Operators::BCs> op_bc_p_;
  Teuchos::RCP<Operators::BCs> op_bc_pc_;
  Teuchos::RCP<Operators::BCs> op_bc_xw_;
  Teuchos::RCP<Operators::BCs> op_bc_xn_;

  // upwind operator
  Teuchos::RCP<Operators::UpwindFlux<RelativePermeability> > upwind_w_;
  Teuchos::RCP<Operators::UpwindFlux<RelativePermeability> > upwind_n_;

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

  static RegisteredPKFactory<Comp_PK> reg_;

};

}  // namespase Flow
}  // namespace Amanzi

#endif

