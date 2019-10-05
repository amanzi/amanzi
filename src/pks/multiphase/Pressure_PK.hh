/*
This is the multiphase flow component of the Amanzi code. 

Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided in the top-level COPYRIGHT file.

Authors: Quan Bui (mquanbui@math.umd.edu)
*/

#ifndef AMANZI_PRESSURE_PK_HH_
#define AMANZI_PRESSURE_PK_HH_

// TPLs
#include "Teuchos_RCP.hpp"

// Basic data structure include
#include "Tensor.hh"
#include "TreeVector.hh"

// Time integration include
#include "PK_PhysicalBDF.hh"

// General include
#include "State.hh"
#include "PDE_DiffusionFactory.hh"
#include "PK_Factory.hh"
#include "PDE_DiffusionFV.hh"

// Specific include for this PK
#include "RelativePermeability.hh"
#include "MultiphaseTypeDefs.hh"
#include "WaterRetentionModel.hh"

#include "LinearOperatorFactory.hh"
#include "FlowBoundaryFunction.hh"
#include "PK_DomainFunction.hh"

namespace Amanzi {
namespace Multiphase {
//class State;

class Pressure_PK: public PK_PhysicalBDF {
public:
  Pressure_PK(Teuchos::ParameterList& pk_tree,
                    const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                    const Teuchos::RCP<State>& S,
                    const Teuchos::RCP<TreeVector>& soln);

  ~Pressure_PK();

  // New interface for a PK
  virtual void Setup();
  virtual void Initialize() {
      InitializeFields();
      InitializePressure();
  }

  virtual double get_dt(){ return dt_; }
  virtual void set_dt(double){};
  virtual bool AdvanceStep(double t_old, double t_new, bool reinit);
  virtual void CommitStep(double t_old, double t_new);
  virtual void CalculateDiagnostics(){};
  virtual std::string name(){return "pressure";}

  // Main methods of this PK
  void InitializeFields();
  void InitializePressure();
  void CommitState(const Teuchos::Ptr<State>& S);

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
  		     Teuchos::RCP<const TreeVector> du) { return 0.0; }

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
                   Teuchos::RCP<TreeVector> du) {
    return AmanziSolvers::FnBaseDefs::CORRECTION_NOT_MODIFIED;
  }

  // experimental approach -- calling this indicates that the time
  // integration scheme is changing the value of the solution in
  // state.
  virtual void ChangedSolution() {}

  // methods to compute boundary and source terms
  void SetAbsolutePermeabilityTensor(Teuchos::RCP<CompositeVector> Sw);
  void AddSourceTerms(CompositeVector& rhs);
  void ComputeBCs();

  // io members, implemented in Pressure_PK_IO.cc
  void ProcessParameterList(Teuchos::ParameterList& list);
  void ProcessStringSourceDistribution(const std::string name, int* method);

  // accesss methods
  Teuchos::RCP<Operators::PDE_Diffusion> op_prec() { return op_preconditioner_; }

public:
  int ncells_owned, ncells_wghost;
  int nfaces_owned, nfaces_wghost;
  int missed_bc_faces_, dirichlet_bc_faces_;

  Teuchos::RCP<Teuchos::ParameterList> linear_operator_list_;
  Teuchos::RCP<Teuchos::ParameterList> preconditioner_list_;
  Teuchos::RCP<Teuchos::ParameterList> mp_list_;

private:
  // mesh structure and geometry
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  int dim_;

  // Stationary physical quantatities
  std::vector<WhetStone::Tensor> K_;
  Teuchos::RCP<Epetra_Vector> Kxy;
  AmanziGeometry::Point gravity_;
  double g_, atm_pressure_, dt_;
  std::vector<double> rho_, mu_;

  bool include_gravity_;

  // source and sink terms
  std::vector<PK_DomainFunction*> srcs_;
  //int src_sink_distribution_;

  // various string objects for names of solver, preconditioners, etc.
  std::string passwd_, solver_name_, preconditioner_name_;

  // Verbose control
  VerboseObject* vo_;

  // boundary conditons
  std::vector<int> bc_submodel;

  Flow::FlowBoundaryFunction* bc_pressure;
  Flow::FlowBoundaryFunction* bc_flux;

  // State and operators
  Teuchos::RCP<State> S_;

  Teuchos::RCP<CompositeVector> tot_mobility_;
  Teuchos::RCP<CompositeVector> gravity_factor_;
  Teuchos::RCP<CompositeVector> rhs_;

  Teuchos::RCP<RelativePermeability> rel_perm_w_;
  Teuchos::RCP<RelativePermeability> rel_perm_n_;
  Teuchos::RCP<Operators::Operator> op_;
  Teuchos::RCP<Operators::PDE_Diffusion> op_matrix_;
  Teuchos::RCP<Operators::PDE_DiffusionFV> op_preconditioner_;
  Teuchos::RCP<Operators::BCs> op_bc_;
  AmanziSolvers::LinearOperatorFactory<Operators::PDE_Diffusion, CompositeVector, CompositeVectorSpace> sfactory;
  Teuchos::RCP<AmanziSolvers::LinearOperator<Operators::PDE_Diffusion, CompositeVector, CompositeVectorSpace> > solver_;

  // The solution obtained from solving for pressure
  Teuchos::RCP<CompositeVector> sol_;
  Teuchos::RCP<TreeVector> sol_tree_;

  static RegisteredPKFactory<Pressure_PK> reg_;

};

}  // namespase Flow
}  // namespace Amanzi

#endif

