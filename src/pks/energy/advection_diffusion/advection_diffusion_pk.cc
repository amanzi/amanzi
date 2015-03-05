/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon
------------------------------------------------------------------------- */

#include "advection_factory.hh"
#include "energy_bc_factory.hh"
#include "CompositeVectorSpace.hh"

#include "OperatorDiffusionMFD.hh"
#include "advection_diffusion.hh"

namespace Amanzi {
namespace Energy {

void AdvectionDiffusion::setup(const Teuchos::Ptr<State>& S) {
  PKPhysicalBDFBase::setup(S);

  // require fields
  Teuchos::RCP<CompositeVectorSpace> factory;

  // -- temperature -- on cells and faces
  std::vector<AmanziMesh::Entity_kind> locations2(2);
  std::vector<std::string> names2(2);
  std::vector<int> num_dofs2(2,1);
  locations2[0] = AmanziMesh::CELL;
  locations2[1] = AmanziMesh::FACE;
  names2[0] = "cell";
  names2[1] = "face";

  factory = S->RequireField("temperature", "energy");
  factory->SetMesh(S->GetMesh());
  factory->SetGhosted(true);
  factory->SetComponents(names2, locations2, num_dofs2);

  // -- thermal conductivity -- just cells
  factory = S->RequireField("thermal_conductivity", "energy");
  factory->SetMesh(S->GetMesh());
  factory->SetGhosted(true);
  factory->SetComponent("cell", AmanziMesh::CELL, 1);

  // -- independent variables (not owned by this pk)
  factory = S->RequireField("darcy_flux");
  factory->SetMesh(S->GetMesh());
  factory->SetGhosted(true);
  factory->SetComponent("face", AmanziMesh::FACE, 1);

  factory = S->RequireField("cell_volume");
  factory->SetMesh(S->GetMesh());
  factory->SetGhosted(true);
  factory->AddComponent("cell", AmanziMesh::CELL, 1);
  S->RequireFieldEvaluator("cell_volume");

  // boundary conditions
  Teuchos::ParameterList bc_plist = plist_->sublist("boundary conditions", true);
  EnergyBCFactory bc_factory(S->GetMesh(), bc_plist);
  bc_temperature_ = bc_factory.CreateTemperature();
  bc_flux_ = bc_factory.CreateEnthalpyFlux();

  int nfaces = S->GetMesh()->num_entities(AmanziMesh::FACE, AmanziMesh::USED);
  bc_markers_.resize(nfaces, Operators::OPERATOR_BC_NONE);
  bc_values_.resize(nfaces, 0.0);
  std::vector<double> mixed;
  bc_ = Teuchos::rcp(new Operators::BCs(Operators::OPERATOR_BC_TYPE_FACE, bc_markers_, bc_values_, mixed));

  // operator for advection terms
  Operators::AdvectionFactory advection_factory;
  Teuchos::ParameterList advect_plist = plist_->sublist("Advection");
  advection_ = advection_factory.create(advect_plist, S->GetMesh());
  advection_->set_num_dofs(1);

  // operator for the diffusion terms
  Teuchos::ParameterList mfd_plist = plist_->sublist("Diffusion");
  matrix_diff_ = Teuchos::rcp(new Operators::OperatorDiffusionMFD(mfd_plist, mesh_));
  matrix_diff_->SetBCs(bc_);
  matrix_diff_->Setup(Teuchos::null);
  matrix_ = matrix_diff_->global_operator();

  // preconditioner
  Teuchos::ParameterList mfd_pc_plist = plist_->sublist("Diffusion PC");
  preconditioner_diff_ = Teuchos::rcp(new Operators::OperatorDiffusionMFD(mfd_pc_plist, mesh_));
  preconditioner_diff_->SetBCs(bc_);
  preconditioner_diff_->Setup(Teuchos::null);
  preconditioner_ = preconditioner_diff_->global_operator();
  preconditioner_acc_ = Teuchos::rcp(new Operators::OperatorAccumulation(AmanziMesh::CELL, preconditioner_));
  preconditioner_->SymbolicAssembleMatrix();
  preconditioner_->InitPreconditioner("preconditioner", mfd_pc_plist);

};


// -- Initialize owned (dependent) variables.
void AdvectionDiffusion::initialize(const Teuchos::Ptr<State>& S) {
  PKPhysicalBDFBase::initialize(S);

  double time = S->time();
  bc_temperature_->Compute(time);
  bc_flux_->Compute(time);
  UpdateBoundaryConditions_();
};


// Evaluate BCs
void AdvectionDiffusion::UpdateBoundaryConditions_() {
  for (unsigned int n=0; n!=bc_markers_.size(); ++n) {
    bc_markers_[n] = Operators::OPERATOR_BC_NONE;
    bc_values_[n] = 0.0;
  }

  Functions::BoundaryFunction::Iterator bc;
  for (bc=bc_temperature_->begin(); bc!=bc_temperature_->end(); ++bc) {
    int f = bc->first;
    bc_markers_[f] = Operators::OPERATOR_BC_DIRICHLET;
    bc_values_[f] = bc->second;
  }

  for (bc=bc_flux_->begin(); bc!=bc_flux_->end(); ++bc) {
    int f = bc->first;
    bc_markers_[f] = Operators::OPERATOR_BC_NEUMANN;
    bc_values_[f] = bc->second;
  }
};


// Add a boundary marker to owned faces.
void AdvectionDiffusion::ApplyBoundaryConditions_(const Teuchos::RCP<CompositeVector>& temperature) {
  Epetra_MultiVector& temp_f = *temperature->ViewComponent("face",true);
  int nfaces = temperature->size("face",false);
  for (int f=0; f!=nfaces; ++f) {
    if (bc_markers_[f] == Operators::OPERATOR_BC_DIRICHLET) {
      temp_f[0][f] = bc_values_[f];
    }
  }
};

} // namespace Energy
} // namespace Amanzi
