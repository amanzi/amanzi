/*
  MultiPhase

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Quan Bui (mquanbui@math.umd.edu)
*/

#ifndef AMANZI_REDUCED2P2C_PK_HH_
#define AMANZI_REDUCED2P2C_PK_HH_

// Amanzi
#include "BDF1_TI.hh"
#include "CombinativeTreeOperator.hh"
#include "LinearOperatorFactory.hh"
#include "TimerManager.hh"
#include "TreeOperator.hh"

// Amanzi::Multiphase
#include "CompW_PK.hh"
#include "CompH_PK.hh"
#include "GasConstraint.hh"

namespace Amanzi {
namespace Multiphase {

class Reduced2p2c_PK: public PK_PhysicalBDF {
public:
  Reduced2p2c_PK(Teuchos::ParameterList& pk_tree,
                 const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                 const Teuchos::RCP<State>& S,
                 const Teuchos::RCP<TreeVector>& soln);

  ~Reduced2p2c_PK() {};

  // New interface for a PK
  virtual void Setup(const Teuchos::Ptr<State>& S) override {};
  virtual void Initialize(const Teuchos::Ptr<State>& S) override;

  virtual double get_dt() override { return dT; }
  virtual void set_dt(double) override {};
  virtual bool AdvanceStep(double t_old, double t_new, bool reinit) override;
  virtual void CommitStep(double t_old, double t_new, const Teuchos::RCP<State>& S) override;
  virtual void CalculateDiagnostics(const Teuchos::RCP<State>& S) override {};
  virtual std::string name() override { return "multiphase multicomponent"; }

  // Time integration interface new_mpc, implemented in Pressure_PK_TI.cc
  // computes the non-linear functional f = f(t,u,udot)
  virtual void FunctionalResidual(double t_old, double t_new, 
                                  Teuchos::RCP<TreeVector> u_old,
                                  Teuchos::RCP<TreeVector> u_new,
                                  Teuchos::RCP<TreeVector> f) override;

  // applies preconditioner to u and returns the result in Pu
  virtual int ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, 
                                  Teuchos::RCP<TreeVector> Pu) override;

  // applies preconditioner to u and returns the result in Pu
  int ApplyJacobian(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu);
  // updates the preconditioner
  virtual void UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h) override;

  // computes a norm on u-du and returns the result
  virtual double ErrorNorm(Teuchos::RCP<const TreeVector> u,
                           Teuchos::RCP<const TreeVector> du) override; 

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
                   Teuchos::RCP<TreeVector> du) override; 

  void ClipSaturation(Teuchos::RCP<CompositeVector> s, double tol);
  void ClipConcentration(Teuchos::RCP<CompositeVector> rho);

  virtual int ReportStatistics() override { return ln_itrs_; }

  // experimental approach -- calling this indicates that the time
  // integration scheme is changing the value of the solution in
  // state.
  virtual void ChangedSolution() override {};

  void ProcessSublistTimeIntegration(Teuchos::ParameterList& list, const std::string name, TI_Specs& ti_specs);
  void ProcessSublistTimeInterval(Teuchos::ParameterList& ti_list,  TI_Specs& ti_specs);
  void ProcessStringTimeIntegration(const std::string name, int* method);
  void ProcessStringErrorOptions(Teuchos::ParameterList& list, int* control);
  std::string FindStringLinearSolver(const Teuchos::ParameterList& plist);
  std::string FindStringPreconditioner(const Teuchos::ParameterList& list);

  // access methods
  Teuchos::RCP<Operators::TreeOperator> op_tree() {return tree_op_;}
  Teuchos::RCP<TreeVector> rhs() {return rhs_;}

public:
  double dT, dT_actual, dTnext, T_physics;
  Teuchos::RCP<Teuchos::ParameterList> ti_list_;
  Teuchos::RCP<Teuchos::ParameterList> linear_operator_list_;
  Teuchos::RCP<Teuchos::ParameterList> pc_list_;
  Teuchos::RCP<Teuchos::ParameterList> cpr_list_;
  Teuchos::ParameterList& pk_tree_;
  Teuchos::RCP<Teuchos::ParameterList> glist_, mpmc_list_;
  double accumulateSolveTime_, accumulateSetupTime_;
  
private:
  // time integration phases
  TI_Specs ti_specs_generic_;
  TI_Specs* ti_specs_;
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_; 
  Teuchos::RCP<State> S_;
  std::string passwd_, jacobian_type_, ncp_type_;
  std::string pc_all_name_, linear_solver_name_;
  std::vector<std::string> pc_block_names_;
  std::vector<int> phase_idx;
  double mu_;

  // solution vectors
  Teuchos::RCP<CompositeVector> p1_;
  Teuchos::RCP<CompositeVector> s1_;
  Teuchos::RCP<CompositeVector> rhl_;
  Teuchos::RCP<TreeVector> p1_tree_;
  Teuchos::RCP<TreeVector> s1_tree_;
  Teuchos::RCP<TreeVector> rhl_tree_;
  Teuchos::RCP<TreeVector> soln_;

  Teuchos::RCP<TreeVector> rhs_;

  int error_control_, num_mat_, ln_itrs_, nl_itrs_, ts_cnt_;

  // verbose object
  VerboseObject* vo_;

  bool cpr_enhanced_;

  Teuchos::RCP<Operators::Operator> op_diff_;
  Teuchos::RCP<Operators::TreeOperator> tree_op_;
  Teuchos::RCP<Operators::TreeOperator> tree_op_precond_;
  Teuchos::RCP<Operators::CombinativeTreeOperator> comb_tree_op_;
  Teuchos::RCP<CompW_PK> comp_w_pk_;
  Teuchos::RCP<CompH_PK> comp_h_pk_;
  Teuchos::RCP<GasConstraint> gas_constraint_pk_;

  Teuchos::RCP<BDF1_TI<TreeVector, TreeVectorSpace> > bdf1_dae;  // Time integrators

  Teuchos::RCP<AmanziSolvers::LinearOperator<Operators::CombinativeTreeOperator, TreeVector, TreeVectorSpace> >
      solver_comb_;
  Teuchos::RCP<AmanziSolvers::LinearOperator<Operators::TreeOperator, TreeVector, TreeVectorSpace> >
      solver_tree_;

  AmanziSolvers::LinearOperatorFactory<Operators::CombinativeTreeOperator, TreeVector, TreeVectorSpace> factory_comb;
  AmanziSolvers::LinearOperatorFactory<Operators::TreeOperator, TreeVector, TreeVectorSpace> factory_tree;

  static RegisteredPKFactory<Reduced2p2c_PK> reg_;
};

}  // namespase Multiphase
}  // namespace Amanzi

#endif

