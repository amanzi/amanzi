/*
  This is the energy component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  This is a derived abstract class.
*/

#ifndef AMANZI_ENERGY_PK_HH_
#define AMANZI_ENERGY_PK_HH_

#include "Epetra_Vector.h"
#include "Epetra_FECrsMatrix.h"
#include "Teuchos_RCP.hpp"

#include "BDFFnBase.hh"
#include "CompositeVectorSpace.hh"
#include "OperatorDiffusionFactory.hh"
#include "VerboseObject.hh"

#include "PK.hh"
#include "primary_variable_field_evaluator.hh"
#include "tensor.hh"

namespace Amanzi {
namespace Energy {

class Energy_PK : public PK, public Amanzi::BDFFnBase<CompositeVector> {
 public:
  Energy_PK(Teuchos::ParameterList& glist, Teuchos::RCP<State> S);
  virtual ~Energy_PK() {};

  // main PK methods
  void Setup() {};
  void Initialize();

  bool AdvanceStep(double t_old, double t_new) {};
  void CommitStep(double t_old, double t_new) {};
  void CalculateDiagnostics() {};

  double get_dt() {};
  void SetState(const Teuchos::RCP<State>& S) { S_ = S; }
  std::string name() { return "energy"; }

  // main energy methods
  void InitializeFields();

  // methods required for time integration (EMPTY so far)
  void Functional(const double Told, double Tnew,
                  Teuchos::RCP<CompositeVector> u_old, Teuchos::RCP<CompositeVector> u_new,
                  Teuchos::RCP<CompositeVector> f) {};
  void ApplyPreconditioner(Teuchos::RCP<const CompositeVector> u, Teuchos::RCP<CompositeVector> Hu) {};
  void UpdatePreconditioner(double T, Teuchos::RCP<const CompositeVector> u, double dT) {};

  double ErrorNorm(Teuchos::RCP<const CompositeVector> u, Teuchos::RCP<const CompositeVector> du) {};
  bool IsAdmissible(Teuchos::RCP<const CompositeVector> up) { return true; }
  bool ModifyPredictor(double dT, Teuchos::RCP<const CompositeVector> u0,
                       Teuchos::RCP<CompositeVector> u) {};
  AmanziSolvers::FnBaseDefs::ModifyCorrectionResult
      ModifyCorrection(double dT, Teuchos::RCP<const CompositeVector> res,
                       Teuchos::RCP<const CompositeVector> u,
                       Teuchos::RCP<CompositeVector> du) {};
  void ChangedSolution() {};

 public:
  int ncells_owned, ncells_wghost;
  int nfaces_owned, nfaces_wghost;

 protected:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  int dim;

  Teuchos::RCP<State> S_;
  std::string passwd_;

 private:
  Teuchos::RCP<Operators::OperatorDiffusion> op_matrix_, op_preconditioner_;
  Teuchos::RCP<Operators::BCs> op_bc_;

  std::vector<WhetStone::Tensor> K; 

  std::vector<int> bc_model_, bc_submodel_; 
  std::vector<double> bc_value_, bc_mixed_; 

 protected:
  Teuchos::ParameterList plist_;
  VerboseObject* vo_;
};

}  // namespace Energy
}  // namespace Amanzi

#endif
