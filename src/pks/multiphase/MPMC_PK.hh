/*
This is the multiphase flow component of the Amanzi code. 

Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided in the top-level COPYRIGHT file.

Authors: Quan Bui (mquanbui@math.umd.edu)
*/

#ifndef AMANZI_MPMC_PK_HH_
#define AMANZI_MPMC_PK_HH_

// Specific include for this PK
#include "Comp_PK.hh"
#include "Phase_Constraint_PK.hh"
#include "TreeOperator.hh"
#include "BDF1_TI.hh"

namespace Amanzi {
namespace Multiphase {
//class State;

class MPMC_PK: public FnTimeIntegratorPK {
public:
  MPMC_PK(Teuchos::ParameterList& pk_tree,
                    const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                    const Teuchos::RCP<State>& S,
                    const Teuchos::RCP<TreeVector>& soln);

  ~MPMC_PK();

  // New interface for a PK
  virtual void Setup(){};
  virtual void Initialize();

  virtual double get_dt(){};
  virtual void set_dt(double){};
  virtual bool AdvanceStep(double t_old, double t_new, bool reinit);
  virtual void CommitStep(double t_old, double t_new);
  virtual void CalculateDiagnostics(){};
  virtual std::string name(){return "multiphase multicomponent";}

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
         Teuchos::RCP<TreeVector> du) {
  }

  // experimental approach -- calling this indicates that the time
  // integration scheme is changing the value of the solution in
  // state.
  virtual void ChangedSolution() {}

  void ProcessSublistTimeIntegration(Teuchos::ParameterList& list, const std::string name, TI_Specs& ti_specs);
  void ProcessSublistTimeInterval(Teuchos::ParameterList& ti_list,  TI_Specs& ti_specs);
  void ProcessStringSourceDistribution(const std::string name, int* method);
  void ProcessStringTimeIntegration(const std::string name, int* method);
  void ProcessStringErrorOptions(Teuchos::ParameterList& list, int* control);
  std::string FindStringLinearSolver(const Teuchos::ParameterList& plist);
  std::string FindStringPreconditioner(const Teuchos::ParameterList& list);

  // access methods
  Teuchos::RCP<Operators::TreeOperator> op_tree() {return tree_op_;}
  Teuchos::RCP<TreeVector> rhs() {return rhs_;}

public:
  double dT, dT_actual, dTnext, T_physics;
  Teuchos::ParameterList ti_list_;
  Teuchos::ParameterList linear_operator_list_;
  Teuchos::ParameterList preconditioner_list_;
  Teuchos::ParameterList& pk_tree_;
  Teuchos::RCP<Teuchos::ParameterList> glist_, mpmc_list_;
  
private:
  // time integration phases
  TI_Specs ti_specs_generic_;
  TI_Specs* ti_specs_;
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_; 
  Teuchos::RCP<State> S_;
  std::string passwd_, jacobian_type_;

  // solution vectors
  Teuchos::RCP<CompositeVector> p1_;
  Teuchos::RCP<CompositeVector> s2_;
  Teuchos::RCP<CompositeVector> fuga1_;
  Teuchos::RCP<CompositeVector> fuga2_;
  Teuchos::RCP<TreeVector> p1_tree_;
  Teuchos::RCP<TreeVector> s2_tree_;
  Teuchos::RCP<TreeVector> fuga1_tree_;
  Teuchos::RCP<TreeVector> fuga2_tree_;
  Teuchos::RCP<TreeVector> soln_;

  Teuchos::RCP<TreeVector> rhs_;

  int error_control_, num_mat_;

  Teuchos::RCP<Operators::Operator> op_diff_;
  //Teuchos::RCP<Operators::TreeOperator> tree_op_;
  //Operators::TreeOperator* tree_op_;
  Teuchos::RCP<Operators::TreeOperator> tree_op_;
  Teuchos::RCP<Comp_PK> comp_w_pk_;
  Teuchos::RCP<Comp_PK> comp_h_pk_;
  Teuchos::RCP<Phase_Constraint_PK> phase_constraint_pk_;

  Teuchos::RCP<BDF1_TI<TreeVector, TreeVectorSpace> > bdf1_dae;  // Time integrators

  static RegisteredPKFactory<MPMC_PK> reg_;

};

}  // namespase Flow
}  // namespace Amanzi

#endif

