/*
This is the multiphase flow component of the Amanzi code. 

Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided in the top-level COPYRIGHT file.

Authors: Quan Bui (mquanbui@math.umd.edu)
*/

#ifndef AMANZI_PHASE1_PK_HH_
#define AMANZI_PHASE1_PK_HH_

// Trilinos include
#include "Teuchos_RCP.hpp"

// Basic data structure include
#include "Tensor.hh"
#include "TreeVector.hh"

// Time integration include
#include "FnTimeIntegratorPK.hh"
//#include "TI_Specs.hh"

// General include
#include "State.hh"
#include "OperatorDiffusionFactory.hh"
#include "OperatorAdvection.hh"
#include "OperatorAccumulation.hh"
#include "Op.hh"
#include "PK_Factory.hh"
#include "primary_variable_field_evaluator.hh"
#include "UpwindFactory.hh"
#include "Upwind.hh"

// Specific include for this PK
#include "RelativePermeability.hh"
#include "MultiphaseTypeDefs.hh"
#include "WaterRetentionModel.hh"

#include "FlowBoundaryFunction.hh"
#include "FlowDomainFunction.hh"
#include "Flow_BC_Factory.hh"
#include "Flow_SourceFactory.hh"

namespace Amanzi {
namespace Multiphase {
//class State;

class Op;

class Phase1_PK: public FnTimeIntegratorPK {
public:
  Phase1_PK(Teuchos::ParameterList& pk_tree,
                    const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                    const Teuchos::RCP<State>& S,
                    const Teuchos::RCP<TreeVector>& soln);

  ~Phase1_PK();

  // New interface for a PK
  virtual void Setup(){};
  virtual void Initialize() {
      InitializeFields();
      InitializePhase1();
      InitNextTI();
  }

  virtual double get_dt(){return 1;}
  virtual void set_dt(double){};
  virtual bool AdvanceStep(double t_old, double t_new, bool reinit);
  virtual void CommitStep(double t_old, double t_new);
  virtual void CalculateDiagnostics(){};
  virtual std::string name(){return "phase1 pk";}

  // Main methods of this PK
  void InitializeFields();
  void InitializePhase1();
  void InitNextTI();

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

  void NumericalJacobian(double t_old, double t_new, 
                         Teuchos::RCP<const TreeVector> u, double eps);

  // methods to compute boundary and source terms
  //void SetAbsolutePermeabilityTensor(Teuchos::RCP<CompositeVector> Sw);
  void SetAbsolutePermeabilityTensor();
  void AddSourceTerms(CompositeVector& rhs) {};
  void ComputeBCs();

  // io members, implemented in Multiphase_IO.cc
  void ProcessParameterList(Teuchos::ParameterList& list);
  void ProcessStringSourceDistribution(const std::string name, int* method);

  // misc function
  void ComputePc(const CompositeVector& Sw, CompositeVector& Pc);
  void AddGradientCapillaryPressure(
    const CompositeVector& Pc, const Epetra_Vector& trans, CompositeVector& rhs);


  // access member functions
  Teuchos::RCP<Operators::OperatorDiffusion> op_prec1() { return op1_preconditioner_; }
  Teuchos::RCP<Operators::OperatorAdvection> op_prec2() { return op2_preconditioner_; }
  std::vector<Teuchos::RCP<Operators::OperatorDiffusion> >& Ops() { return ops_; }

public:
  int ncells_owned_, ncells_wghost_;
  int nfaces_owned_, nfaces_wghost_;
  int missed_bc_faces_, dirichlet_bc_faces_;

  Teuchos::RCP<Teuchos::ParameterList> mp_list_;

private:
  // mesh structure and geometry
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  int dim_;

  // Stationary physical quantatities
  std::vector<WhetStone::Tensor> K_;
  Teuchos::RCP<Epetra_Vector> Kxy_;
  AmanziGeometry::Point gravity_;
  double g_, atm_pressure_;
  double rho_, mu_;
  double phi_;

  // source and sink terms
  Flow::FlowDomainFunction* src_sink_;
  int src_sink_distribution_;
  std::string passwd_;

  // Verbose control
  VerboseObject* vo_;

  // boundary conditons
  std::vector<int> bc_model_p_, bc_submodel_p_;
  std::vector<double> bc_value_p_, bc_coef_p_;
  std::vector<int> bc_model_s_, bc_submodel_s_;
  std::vector<double> bc_value_s_, bc_coef_s_;
  std::vector<double> bc_mixed_;

  Flow::FlowBoundaryFunction* bc_pressure_;
  Flow::FlowBoundaryFunction* bc_saturation_;
  Flow::FlowBoundaryFunction* bc_flux_phase1_;
  
  // State
  Teuchos::RCP<State> S_;

  Teuchos::RCP<RelativePermeability> rel_perm_w_;
  Teuchos::RCP<Operators::OperatorDiffusion> op1_matrix_;
  Teuchos::RCP<Operators::OperatorDiffusion> op2_matrix_;
  Teuchos::RCP<Operators::OperatorDiffusion> op1_preconditioner_;
  Teuchos::RCP<Operators::OperatorDiffusion> op_test_;
  Teuchos::RCP<Operators::OperatorAdvection> op2_preconditioner_;
  Teuchos::RCP<Operators::OperatorAccumulation> op_acc_;
  Teuchos::RCP<Operators::BCs> op_bc_s_;
  Teuchos::RCP<Operators::BCs> op_bc_p_;

  // upwind operator
  Teuchos::RCP<Operators::Upwind<RelativePermeability> > upwind_;

  // The solution obtained from solving for pressure
  Teuchos::RCP<CompositeVector> saturation_phase1_;

  std::vector<Teuchos::RCP<Operators::OperatorDiffusion> > ops_;
  typedef std::vector<Teuchos::RCP<Operators::OperatorDiffusion> >::iterator op_iter;
  typedef std::vector<Teuchos::RCP<Operators::Op> >::iterator local_op_iter;

  // solution tree vector
  Teuchos::RCP<TreeVector> soln_;

  // upwind flux
  Teuchos::RCP<CompositeVector> upwind_velocity_;
  Teuchos::RCP<CompositeVector> tmp_flux_; 

  static RegisteredPKFactory<Phase1_PK> reg_;

};

}  // namespase Flow
}  // namespace Amanzi

#endif

