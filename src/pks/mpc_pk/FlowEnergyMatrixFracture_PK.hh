/*
  MPC PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov
           Daniil Svyatskiy

  Process kernel that couples flow and energy in matrix and fractures.
*/

#ifndef AMANZI_FLOW_ENERGY_MATRIX_FRACTURE_PK_HH_
#define AMANZI_FLOW_ENERGY_MATRIX_FRACTURE_PK_HH_

#include "Teuchos_RCP.hpp"

#include "independent_variable_field_evaluator_fromfunction.hh"
#include "secondary_variable_field_evaluator.hh"
#include "PK_BDF.hh"
#include "PK_MPCStrong.hh"
#include "PK_Factory.hh"
#include "TreeOperator.hh"

namespace Amanzi {

class FlowEnergyMatrixFracture_PK : public PK_MPCStrong<PK_BDF> {
 public:
  FlowEnergyMatrixFracture_PK(Teuchos::ParameterList& pk_tree,
                              const Teuchos::RCP<Teuchos::ParameterList>& glist,
                              const Teuchos::RCP<State>& S,
                              const Teuchos::RCP<TreeVector>& soln);

  // PK methods
  virtual void Setup(const Teuchos::Ptr<State>& S);
  virtual void Initialize(const Teuchos::Ptr<State>& S);

  // -- dt is the minimum of the sub pks
  // virtual double get_dt();
  // virtual void set_dt(double dt);

  // -- advance each sub pk from t_old to t_new.
  virtual bool AdvanceStep(double t_old, double t_new, bool reinit = false);
  // virtual void CommitStep(double t_old, double t_new);

  virtual void FunctionalResidual(
      double t_old, double t_new,
      Teuchos::RCP<TreeVector> u_old, Teuchos::RCP<TreeVector> u_new,
      Teuchos::RCP<TreeVector> f);
  
  // updates the preconditioner
  virtual void UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double h);
  
  // // preconditioner application
  virtual int ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu);

  std::string name() { return "thermal flow matrix fracture"; } 

  // virtual void CalculateDiagnostics() {};
  Teuchos::RCP<const Teuchos::ParameterList> linear_operator_list_;
  Teuchos::RCP<const Teuchos::ParameterList> preconditioner_list_;
  Teuchos::RCP<Teuchos::ParameterList> ti_list_;
  
 private:
  const Teuchos::RCP<Teuchos::ParameterList>& glist_;
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_domain_, mesh_fracture_;

  // factory registration
  static RegisteredPKFactory<FlowEnergyMatrixFracture_PK> reg_;
};

}  // namespace Amanzi
#endif
