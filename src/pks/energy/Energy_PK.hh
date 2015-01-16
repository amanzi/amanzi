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

#include "BDFFnBase.hh"
#include "CompositeVector.hh"
#include "OperatorDiffusionFactory.hh"
#include "PK.hh"
#include "primary_variable_field_evaluator.hh"
#include "tensor.hh"
#include "TreeVector.hh"
#include "VerboseObject.hh"

#include "EnergyBoundaryFunction.hh"

namespace Amanzi {
namespace Energy {

class Energy_PK : public PK, public Amanzi::BDFFnBase<TreeVector> {
 public:
  Energy_PK(Teuchos::RCP<const Teuchos::ParameterList>& glist, Teuchos::RCP<State> S);
  virtual ~Energy_PK() {};

  // main PK methods
  virtual void Setup() = 0;
  virtual void Initialize();

  bool AdvanceStep(double t_old, double t_new) { return true; }
  void CommitStep(double t_old, double t_new) {};
  void CalculateDiagnostics() {};

  double get_dt() { return 0.0; }
  void SetState(const Teuchos::RCP<State>& S) { S_ = S; }
  std::string name() { return "energy"; }

  // main energy methods
  void InitializeFields();

  // methods required for time integration (EMPTY so far)
  void Functional(const double Told, double Tnew,
                  Teuchos::RCP<TreeVector> u_old, Teuchos::RCP<TreeVector> u_new,
                  Teuchos::RCP<TreeVector> f);
  void ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Hu);
  void UpdatePreconditioner(double T, Teuchos::RCP<const TreeVector> up, double dT);

  double ErrorNorm(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<const TreeVector> du);
  bool IsAdmissible(Teuchos::RCP<const TreeVector> up) {
   return true;
  }
  bool ModifyPredictor(double dT, Teuchos::RCP<const TreeVector> u0, Teuchos::RCP<TreeVector> u) {
    return false;
  }
  AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
      ModifyCorrection(double dT, Teuchos::RCP<const TreeVector> res,
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
  Teuchos::RCP<PrimaryVariableFieldEvaluator> temperature_eval;

  Teuchos::RCP<State> S_;
  std::string passwd_;

  // conductivity tensor
  std::vector<WhetStone::Tensor> K; 

  // boundary conditons
  EnergyBoundaryFunction* bc_temperature; 
  EnergyBoundaryFunction* bc_flux; 

  std::vector<int> bc_model_, bc_submodel_; 
  std::vector<double> bc_value_, bc_mixed_; 
  int dirichlet_bc_faces_;

  // operators
  Teuchos::RCP<Operators::OperatorDiffusion> op_matrix_, op_preconditioner_;
  Teuchos::RCP<Operators::BCs> op_bc_;

 protected:
  Teuchos::RCP<const Teuchos::ParameterList> glist_;
  VerboseObject* vo_;

  Key energy_key_;  // keys
  Key enthalpy_key_;
  Key conductivity_key_;
  Key uw_conductivity_key_;
};

}  // namespace Energy
}  // namespace Amanzi

#endif
