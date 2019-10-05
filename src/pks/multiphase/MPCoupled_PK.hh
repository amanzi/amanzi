/*
This is the multiphase flow component of the Amanzi code. 

Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided in the top-level COPYRIGHT file.

Authors: Quan Bui (mquanbui@math.umd.edu)
*/

#ifndef AMANZI_MPCOUPLED_PK_HH_
#define AMANZI_MPCOUPLED_PK_HH_

// Specific include for this PK
#include "Phase1_PK.hh"
#include "Phase2_PK.hh"
#include "TotalPhase_PK.hh"
#include "CombinativeTreeOperator.hh"
#include "BDF1_TI.hh"
#include "TimerManager.hh"

namespace Amanzi {
namespace Multiphase {
//class State;

class MPCoupled_PK: public FnTimeIntegratorPK {
public:
  MPCoupled_PK(Teuchos::ParameterList& pk_tree,
                    const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                    const Teuchos::RCP<State>& S,
                    const Teuchos::RCP<TreeVector>& soln);

  ~MPCoupled_PK();

  // New interface for a PK
  virtual void Setup(){};
  virtual void Initialize();

  virtual double get_dt(){ return dt_; };
  virtual void set_dt(double dt){ dt_ = dt; };
  virtual bool AdvanceStep(double t_old, double t_new, bool reinit);
  virtual void CommitStep(double t_old, double t_new);
  virtual void CalculateDiagnostics(){};
  virtual std::string name(){return "multiphase coupled";}

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
           Teuchos::RCP<const TreeVector> du); 

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
         Teuchos::RCP<TreeVector> du);

  virtual int ReportStatistics() { return ln_itrs_; }

  // experimental approach -- calling this indicates that the time
  // integration scheme is changing the value of the solution in
  // state.
  virtual void ChangedSolution() {}

  void ClipSaturation(Teuchos::RCP<CompositeVector> s, double tol);

  void ComputeAPinv();

  // access methods
  Teuchos::RCP<const Operators::TreeOperator> op_tree() {return tree_op_;}

public:
  double dT_actual, dTnext, T_physics;
  Teuchos::RCP<Teuchos::ParameterList> ti_list_;
  Teuchos::RCP<Teuchos::ParameterList> linear_operator_list_;
  Teuchos::RCP<Teuchos::ParameterList> pc_list_;
  Teuchos::RCP<Teuchos::ParameterList> cpr_list_;
  Teuchos::ParameterList& pk_tree_;
  Teuchos::RCP<Teuchos::ParameterList> glist_;
  
private:
  // Epetra communicator
  Epetra_MpiComm* comm_;
  int *coarse_indices_array_;

  // time integration phases
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_; 
  Teuchos::RCP<State> S_;
  std::string passwd_, jacobian_type_;
  std::string pc_all_name_, pc_ass_name_, linear_solver_name_;
  std::vector<std::string> pc_block_names_;
  bool cpr_enhanced_, block_factorization_;
  bool include_capillary_, include_gravity_, use_total_mobility_;

  double block10_scaling_, block01_scaling_;

  // internal time step
  double dt_;

  // verbose object
  VerboseObject* vo_;
  
  // solution vectors
  Teuchos::RCP<CompositeVector> p1_;
  Teuchos::RCP<CompositeVector> s2_;
  Teuchos::RCP<CompositeVector> s1_;
  Teuchos::RCP<TreeVector> p1_tree_;
  Teuchos::RCP<TreeVector> s2_tree_;
  Teuchos::RCP<TreeVector> soln_;

  int nl_itrs_, ln_itrs_, max_ln_itrs_, min_ln_itrs_, num_mat_;

  Teuchos::RCP<Operators::Operator> op_diff_;
  Teuchos::RCP<Operators::TreeOperator> tree_op_;
  Teuchos::RCP<Operators::CombinativeTreeOperator> comb_tree_op_;
  Teuchos::RCP<Phase1_PK> phase1_pk_;
  Teuchos::RCP<Phase2_PK> phase2_pk_;
  Teuchos::RCP<TotalPhase_PK> tot_phase_pk_;

  Teuchos::RCP<BDF1_TI<TreeVector, TreeVectorSpace> > bdf1_dae;  // Time integrators
  Teuchos::RCP<AmanziSolvers::LinearOperator<Operators::CombinativeTreeOperator, TreeVector, TreeVectorSpace> >
      solver_comb_;
  Teuchos::RCP<AmanziSolvers::LinearOperator<Operators::TreeOperator, TreeVector, TreeVectorSpace> >
      solver_tree_;

  AmanziSolvers::LinearOperatorFactory<Operators::CombinativeTreeOperator, TreeVector, TreeVectorSpace> factory_comb;
  AmanziSolvers::LinearOperatorFactory<Operators::TreeOperator, TreeVector, TreeVectorSpace> factory_tree;

  static RegisteredPKFactory<MPCoupled_PK> reg_;
};

}  // namespase Flow
}  // namespace Amanzi

#endif

