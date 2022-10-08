/*
  MPC PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov
           Daniil Svyatskiy

  Process kernel that couples Darcy flow in matrix and fracture network.
*/

#ifndef AMANZI_DARCY_MATRIX_FRACTURE_PK_HH_
#define AMANZI_DARCY_MATRIX_FRACTURE_PK_HH_

#include "Teuchos_RCP.hpp"

#include "EvaluatorIndependentFunction.hh"
#include "EvaluatorSecondary.hh"
#include "PK_BDF.hh"
#include "PK_MPCStrong.hh"
#include "PK_Factory.hh"
#include "TreeOperator.hh"

namespace Amanzi {

class FlowMatrixFracture_PK : public PK_MPCStrong<PK_BDF> {
 public:
  FlowMatrixFracture_PK(Teuchos::ParameterList& pk_tree,
                        const Teuchos::RCP<Teuchos::ParameterList>& glist,
                        const Teuchos::RCP<State>& S,
                        const Teuchos::RCP<TreeVector>& soln);

  // PK methods
  virtual void Setup() override;
  virtual void Initialize() override;

  // -- dt is the minimum of the sub pks
  // virtual double get_dt();
  // virtual void set_dt(double dt);

  // -- advance each sub pk from t_old to t_new.
  virtual bool AdvanceStep(double t_old, double t_new, bool reinit = false) override;

  virtual
  void FunctionalResidual(double t_old, double t_new, Teuchos::RCP<TreeVector> u_old,
                          Teuchos::RCP<TreeVector> u_new, Teuchos::RCP<TreeVector> f) override;
  
  // updates the preconditioner
  virtual void UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> u, double dt) override;
  
  // // preconditioner application
  virtual
  int ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu) override;

  std::string name() override { return "flow matrix-fracture"; } 

  // virtual void CalculateDiagnostics() {};
  Teuchos::RCP<const Teuchos::ParameterList> linear_operator_list_;
  Teuchos::RCP<const Teuchos::ParameterList> preconditioner_list_;
  Teuchos::RCP<Teuchos::ParameterList> ti_list_;
  
 private:
  void AddDefaultPrimaryEvaluator_(const Key& key, const Tag& tag);

 private:
  const Teuchos::RCP<Teuchos::ParameterList>& glist_;
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_domain_, mesh_fracture_;

  Teuchos::RCP<Operators::TreeOperator> op_matrix_, op_preconditioner_;

  // factory registration
  static RegisteredPKFactory<FlowMatrixFracture_PK> reg_;
};

}  // namespace Amanzi
#endif
