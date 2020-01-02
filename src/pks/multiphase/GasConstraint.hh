/*
  Multiphase

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Quan Bui (mquanbui@math.umd.edu)
*/

#ifndef AMANZI_GAS_CONSTRAINT_PK_HH_
#define AMANZI_GAS_CONSTRAINT_PK_HH_

// TPLs
#include "Teuchos_RCP.hpp"

// Amanzi
#include "PDE_AdvectionUpwind.hh"
#include "PDE_Accumulation.hh"
#include "PK_Factory.hh"
#include "State.hh"
#include "TreeVector.hh"

// Amanzi::Multiphase
#include "CapillaryPressure.hh"
#include "MultiphaseTypeDefs.hh"

namespace Amanzi {
namespace Multiphase {

class GasConstraint {
 public:
  GasConstraint(Teuchos::ParameterList plist, const Teuchos::RCP<State> S);
  ~GasConstraint() {};

  // New interface for a PK
  void Initialize();
  std::string name() { return "phase constraint pk"; }

  // Main methods of this PK

  // Time integration interface new_mpc, implemented in Pressure_PK_TI.cc
  // computes the non-linear functional f = f(t,u,udot)
  virtual void FunctionalResidual(double t_old, double t_new, 
                                  Teuchos::RCP<TreeVector> u_old,
                                  Teuchos::RCP<TreeVector> u_new,
                                  Teuchos::RCP<TreeVector> f);

  // updates the preconditioner
  virtual void UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h);


  // access methods
  Teuchos::RCP<Operators::PDE_Accumulation> op_prec1() { return op1_acc_; }
  Teuchos::RCP<Operators::PDE_Accumulation> op_prec2() { return op2_acc_; }
  Teuchos::RCP<Operators::PDE_Accumulation> op_prec3() { return op3_acc_; }
  Teuchos::RCP<Operators::PDE_Accumulation> op_prec1_tmp() { return op1_acc_tmp_; }
  Teuchos::RCP<Operators::PDE_Accumulation> op_prec2_tmp() { return op2_acc_tmp_; }
  Teuchos::RCP<Operators::PDE_Accumulation> op_prec3_tmp() { return op3_acc_tmp_; }
  int* getInactiveGasIndices() { return inactive_gas_idx_; }
  int getNumInactiveCells() { return cnt_; }
  void SetNCPFunctionType(std::string ncp_type) { ncp_type_ = ncp_type; }
  void setMu(double mu) { mu_ = mu; }

public:
  Teuchos::ParameterList plist_;

private:
  // mesh structure and geometry
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  int dim_;
  double H_, M_h_, mu_;
  std::string passwd_, ncp_type_;
  int ncells_owned_, ncells_wghost_;
  int nfaces_owned_, nfaces_wghost_;

  // boundary conditons
  std::vector<int> bc_model_;
  std::vector<double> bc_value_;
  std::vector<double> bc_mixed_;
  int* inactive_gas_idx_;
  int cnt_;

  // State and operators
  Teuchos::RCP<State> S_;

  Teuchos::RCP<CapillaryPressure> capillary_pressure_;

  Teuchos::RCP<Operators::PDE_Accumulation> op1_acc_;
  Teuchos::RCP<Operators::PDE_Accumulation> op2_acc_;
  Teuchos::RCP<Operators::PDE_Accumulation> op3_acc_;

  Teuchos::RCP<Operators::PDE_Accumulation> op1_acc_tmp_;
  Teuchos::RCP<Operators::PDE_Accumulation> op2_acc_tmp_;
  Teuchos::RCP<Operators::PDE_Accumulation> op3_acc_tmp_;
  Teuchos::RCP<Operators::BCs> op_bc_;

  // The solution obtained from solving for pressure
  Teuchos::RCP<CompositeVector> sol_;

  // solution tree vector
  Teuchos::RCP<TreeVector> soln_;

  //static RegisteredPKFactory<GasConstraint> reg_;

};

}  // namespase Flow
}  // namespace Amanzi

#endif

