/*
  This is the energy component of the Amanzi code. 
  This is a base class for energy equations.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_ENERGY_PK_HH_
#define AMANZI_ENERGY_PK_HH_

#include "Epetra_Vector.h"
#include "Epetra_FECrsMatrix.h"
#include "Teuchos_RCP.hpp"

#include "CompositeVector.hh"
#include "FnTimeIntegratorPK.hh"
#include "PK.hh"
#include "primary_variable_field_evaluator.hh"
#include "Operator.hh"
#include "OperatorDiffusion.hh"
#include "OperatorAccumulation.hh"
#include "tensor.hh"
#include "TreeVector.hh"
#include "VerboseObject.hh"

#include "EnergyBoundaryFunction.hh"

namespace Amanzi {
namespace Energy {

class Energy_PK : public FnTimeIntegratorPK {
 public:
  Energy_PK(const Teuchos::RCP<Teuchos::ParameterList>& glist, Teuchos::RCP<State> S);
  virtual ~Energy_PK() {};

  // main PK methods
  virtual void Setup();
  virtual void Initialize();

  virtual bool AdvanceStep(double t_old, double t_new) { return true; }
  virtual void CommitStep(double t_old, double t_new) {};
  void CalculateDiagnostics() {};

  double get_dt() { return 0.0; }
  void set_dt(double dt) {}
  void SetState(const Teuchos::RCP<State>& S) { S_ = S; }
  virtual std::string name() { return "energy"; }

  // main energy methods
  void InitializeFields();

  // methods required for time integration
  virtual void ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> hu) {
    op_preconditioner_->ApplyInverse(*u->Data(), *hu->Data());
  }
  bool IsAdmissible(Teuchos::RCP<const TreeVector> up) {
    return true;
  }
  bool ModifyPredictor(double dt, Teuchos::RCP<const TreeVector> u0, Teuchos::RCP<TreeVector> u) {
    return false;
  }
  AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
      ModifyCorrection(double dt, Teuchos::RCP<const TreeVector> res,
                       Teuchos::RCP<const TreeVector> u,
                       Teuchos::RCP<TreeVector> du) {};
  void ChangedSolution() {};

  // other methods
  bool UpdateConductivityData(const Teuchos::Ptr<State>& S);
  void UpdateSourceBoundaryData(double T0, double T1, const CompositeVector& u);
  void ComputeBCs(const CompositeVector& u);

  // access methods for unit tests
  std::vector<WhetStone::Tensor>& get_K() { return K; } 
  Teuchos::RCP<PrimaryVariableFieldEvaluator>& get_temperature_eval() { return temperature_eval; }

 public:
  int ncells_owned, ncells_wghost;
  int nfaces_owned, nfaces_wghost;

 protected:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  int dim;

  const Teuchos::RCP<Teuchos::ParameterList> glist_;
  Teuchos::RCP<Teuchos::ParameterList> ep_list_;
  Teuchos::RCP<const Teuchos::ParameterList> preconditioner_list_;
  Teuchos::RCP<Teuchos::ParameterList> ti_list_;

  // state and primary field
  Teuchos::RCP<State> S_;
  std::string passwd_;
  Teuchos::RCP<PrimaryVariableFieldEvaluator> temperature_eval;

  // keys
  Key energy_key_, prev_energy_key_;
  Key enthalpy_key_;
  Key conductivity_key_;

  // conductivity tensor
  std::vector<WhetStone::Tensor> K; 

  // boundary conditons
  EnergyBoundaryFunction* bc_temperature; 
  EnergyBoundaryFunction* bc_flux; 

  std::vector<int> bc_model_, bc_submodel_; 
  std::vector<double> bc_value_, bc_mixed_; 
  int dirichlet_bc_faces_;

  // operators and solvers
  Teuchos::RCP<Operators::OperatorDiffusion> op_matrix_diff_, op_preconditioner_diff_;
  Teuchos::RCP<Operators::OperatorAccumulation> op_acc_;
  Teuchos::RCP<Operators::Operator> op_matrix_, op_preconditioner_;
  Teuchos::RCP<Operators::BCs> op_bc_;
  std::string preconditioner_name_;

 protected:
  VerboseObject* vo_;
};

}  // namespace Energy
}  // namespace Amanzi

#endif
