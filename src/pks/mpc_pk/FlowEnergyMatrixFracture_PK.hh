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

#include "EvaluatorIndependentFunction.hh"
#include "EvaluatorSecondary.hh"
#include "PDE_CouplingFlux.hh"
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
  virtual void UpdatePreconditioner(double t, Teuchos::RCP<const TreeVector> up, double dt);
  
  // preconditioner application
  virtual int ApplyPreconditioner(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu);

  virtual double ErrorNorm(Teuchos::RCP<const TreeVector> u, 
                            Teuchos::RCP<const TreeVector> du);

  std::string name() { return "thermal flow matrix fracture"; } 

 private:
  // use flag to avoid double counting of coupling terms for Darcy PK
  std::vector<Teuchos::RCP<Operators::PDE_CouplingFlux> > AddCouplingFluxes_(
      Teuchos::RCP<CompositeVectorSpace>& cvs_matrix,
      Teuchos::RCP<CompositeVectorSpace>& cvs_fracture,
      std::shared_ptr<const std::vector<std::vector<int> > > inds_matrix,
      std::shared_ptr<const std::vector<std::vector<int> > > inds_fracture,
      std::shared_ptr<const std::vector<double> > values,
      int i, int j, Teuchos::RCP<Operators::TreeOperator>& op_tree);

  void UpdateCouplingFluxes_(
      const std::vector<Teuchos::RCP<Operators::PDE_CouplingFlux> >& adv_coupling);

  void SwapEvaluatorField_(
      const Key& key, const std::string& passwd,
      Teuchos::RCP<CompositeVector>& fdm_copy,
      Teuchos::RCP<CompositeVector>& fdf_copy);

 public:
  // virtual void CalculateDiagnostics() {};
  Teuchos::RCP<const Teuchos::ParameterList> linear_operator_list_;
  Teuchos::RCP<const Teuchos::ParameterList> preconditioner_list_;
  Teuchos::RCP<Teuchos::ParameterList> ti_list_;
  
 private:
  const Teuchos::RCP<Teuchos::ParameterList>& glist_;
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_domain_, mesh_fracture_;

  Key normal_permeability_key_, normal_conductivity_key_;

  std::vector<Teuchos::RCP<Operators::PDE_CouplingFlux> > adv_coupling_matrix_, adv_coupling_pc_;

  Teuchos::RCP<Matrix<TreeVector, TreeVectorSpace>> op_pc_solver_;

  double residual_norm_;

  // factory registration
  static RegisteredPKFactory<FlowEnergyMatrixFracture_PK> reg_;
};


// non-member function
int ApplyFlattened(const Operators::TreeOperator& op, const TreeVector& X, TreeVector& Y);


// supporting class that is used to due to ApplyFlattened
namespace Operators {

class FlatTreeOperator : public Operators::TreeOperator {
 public:
  using Vector_t = TreeVector;
  using VectorSpace_t = TreeVector::VectorSpace_t;

  FlatTreeOperator(const Teuchos::RCP<const TreeVectorSpace>& tvs) 
    : TreeOperator(tvs) {};

  virtual int Apply(const TreeVector& X, TreeVector& Y) const override {
    return ApplyFlattened(*this, X, Y);
  }
};

}  // namespace Operators

}  // namespace Amanzi
#endif
