/*
This is the multiphase flow component of the Amanzi code. 

Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided in the top-level COPYRIGHT file.

Authors: Quan Bui (mquanbui@math.umd.edu)
*/

#ifndef AMANZI_PHASE_CONSTRAINT_PK_HH_
#define AMANZI_PHASE_CONSTRAINT_PK_HH_

// Trilinos include
#include "Teuchos_RCP.hpp"

// Basic data structure include
#include "TreeVector.hh"

// General include
#include "State.hh"
#include "OperatorAdvection.hh"
#include "OperatorAccumulation.hh"
#include "PK_Factory.hh"

// Specific include for this PK
#include "RelativePermeability.hh"
#include "CapillaryPressure.hh"
#include "MultiphaseTypeDefs.hh"
#include "WaterRetentionModel.hh"

//#include "FlowBoundaryFunction.hh"
//#include "FlowDomainFunction.hh"
//#include "Flow_BC_Factory.hh"
//#include "Flow_SourceFactory.hh"

namespace Amanzi {
namespace Multiphase {

class Phase_Constraint_PK//: public FnTimeIntegratorPK {
{
public:
  Phase_Constraint_PK(Teuchos::ParameterList plist,
                      const Teuchos::RCP<State> S);

  ~Phase_Constraint_PK() {};

  // New interface for a PK
  virtual void Setup(){};
  virtual void Initialize(); 
  virtual std::string name(){return "phase constraint pk";}

  // Main methods of this PK
  int Phase_ID() {return phase_id_;}
  //void CommitState(const Teuchos::Ptr<State>& S);

  // Time integration interface new_mpc, implemented in Pressure_PK_TI.cc
  // computes the non-linear functional f = f(t,u,udot)
  virtual void Functional(double t_old, double t_new, 
                          Teuchos::RCP<TreeVector> u_old,
                          Teuchos::RCP<TreeVector> u_new,
                          Teuchos::RCP<TreeVector> f);

  // updates the preconditioner
  virtual void UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up,
                    double h);


  // access methods
  Teuchos::RCP<Operators::OperatorAdvection> op_prec1() { return op1_prec_; }
  Teuchos::RCP<Operators::OperatorAdvection> op_prec2() { return op2_prec_; }
  Teuchos::RCP<Operators::OperatorAdvection> op_prec3() { return op3_prec_; }
  Teuchos::RCP<Operators::OperatorAdvection> op_prec4() { return op4_prec_; }
  Teuchos::RCP<Operators::OperatorAdvection> op_prec5() { return op5_prec_; }
  Teuchos::RCP<Operators::OperatorAdvection> op_p2_sat_prec() { return op_p2_sat_prec_; }
  Teuchos::RCP<Operators::OperatorAdvection> op_p1_sat_prec() { return op_p1_sat_prec_; }

public:
  Teuchos::ParameterList plist_;

private:
  // mesh structure and geometry
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  int dim_, phase_id_;
  double henry_coef_, P_vap_;
  std::string passwd_;
  int ncells_owned_, ncells_wghost_;
  int nfaces_owned_, nfaces_wghost_;

  // boundary conditons
  std::vector<int> bc_model_, bc_submodel_;
  std::vector<double> bc_value_, bc_coef_;
  std::vector<double> bc_mixed_;

  // State and operators
  Teuchos::RCP<State> S_;

  Teuchos::RCP<RelativePermeability> rel_perm_w_;
  Teuchos::RCP<RelativePermeability> rel_perm_n_;
  Teuchos::RCP<CapillaryPressure> capillary_pressure_;

  Teuchos::RCP<Operators::OperatorAccumulation> op1_acc_;
  Teuchos::RCP<Operators::OperatorAccumulation> op2_acc_;
  Teuchos::RCP<Operators::OperatorAccumulation> op3_acc_;
  Teuchos::RCP<Operators::OperatorAccumulation> op4_acc_;
  Teuchos::RCP<Operators::OperatorAccumulation> op5_acc_;
  Teuchos::RCP<Operators::OperatorAccumulation> op6_acc_;
  Teuchos::RCP<Operators::OperatorAccumulation> op_p1_sat_acc_;
  Teuchos::RCP<Operators::OperatorAccumulation> op_p2_sat_acc_;
  Teuchos::RCP<Operators::OperatorAdvection> op1_prec_;
  Teuchos::RCP<Operators::OperatorAdvection> op2_prec_;
  Teuchos::RCP<Operators::OperatorAdvection> op3_prec_;
  Teuchos::RCP<Operators::OperatorAdvection> op4_prec_;
  Teuchos::RCP<Operators::OperatorAdvection> op5_prec_;
  Teuchos::RCP<Operators::OperatorAdvection> op6_prec_;
  Teuchos::RCP<Operators::OperatorAdvection> op_p1_sat_prec_;
  Teuchos::RCP<Operators::OperatorAdvection> op_p2_sat_prec_;
  Teuchos::RCP<Operators::BCs> op_bc_;

  // The solution obtained from solving for pressure
  Teuchos::RCP<CompositeVector> sol_;

  // solution tree vector
  Teuchos::RCP<TreeVector> soln_;

  static RegisteredPKFactory<Phase_Constraint_PK> reg_;

};

}  // namespase Flow
}  // namespace Amanzi

#endif

