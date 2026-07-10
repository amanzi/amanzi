/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov
           Daniil Svyatskiy
*/

/*
  MPC PK

  Process kernel that couples flow and energy in matrix and fractures.
*/

#ifndef AMANZI_FLOW_ENERGYPH_MATRIX_FRACTURE_PK_HH_
#define AMANZI_FLOW_ENERGYPH_MATRIX_FRACTURE_PK_HH_

#include "Teuchos_RCP.hpp"

#include "EvaluatorIndependentFunction.hh"
#include "EvaluatorSecondary.hh"
#include "PDE_CouplingFlux.hh"
#include "PK_BDF.hh"
#include "PK_MPCStrong.hh"
#include "PK_Factory.hh"
#include "TreeOperator.hh"

namespace Amanzi {

class FlowEnergyPHMatrixFracture_PK : public PK_MPCStrong<PK_BDF> {
 public:
  FlowEnergyPHMatrixFracture_PK(Teuchos::ParameterList& pk_tree,
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

  virtual void FunctionalResidual(double t_old,
                                  double t_new,
                                  Teuchos::RCP<const TreeVector> u_old,
                                  Teuchos::RCP<TreeVector> u_new,
                                  Teuchos::RCP<TreeVector> f) override;

  // -- preconditioner
  virtual void UpdatePreconditioner(double t,
                                    Teuchos::RCP<const TreeVector> up,
                                    double dt) override;

  virtual int ApplyPreconditioner(Teuchos::RCP<const TreeVector> u,
                                  Teuchos::RCP<TreeVector> Pu) override;

  // -- error norm for coupled system
  virtual double ErrorNorm(Teuchos::RCP<const TreeVector> u,
                           Teuchos::RCP<const TreeVector> du) override;

  std::string name() override { return "flow and energy matrix fracture"; }

 public:
  // virtual void CalculateDiagnostics() {};
  Teuchos::RCP<const Teuchos::ParameterList> linear_operator_list_;
  Teuchos::RCP<const Teuchos::ParameterList> preconditioner_list_;
  Teuchos::RCP<Teuchos::ParameterList> ti_list_;

 private:

  void InitializeFlowCoupling_();
  void InitializeHeatCoupling_();
  
  const Teuchos::RCP<Teuchos::ParameterList>& glist_;
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_matrix_, mesh_fracture_;

  Teuchos::RCP<Matrix<TreeVector, TreeVectorSpace>> op_pc_solver_;

  std::vector<Teuchos::RCP<Operators::PDE_CouplingFlux>> thermo_coupling_pc_ops_;
  std::vector<Teuchos::RCP<Operators::PDE_CouplingFlux>> thermo_coupling_mat_ops_;
  
  Key heat_diffusion_to_matrix_key_;
  Key diffusion_to_matrix_key_;
  
  double residual_norm_;

  // factory registration
  static RegisteredPKFactory<FlowEnergyPHMatrixFracture_PK> reg_;
};

} // namespace Amanzi
#endif
