/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon
------------------------------------------------------------------------- */

#include "bdf2_time_integrator.hh"
#include "advection_factory.hh"
#include "energy_bc_factory.hh"
#include "composite_vector_factory.hh"

#include "advection_diffusion.hh"

namespace Amanzi {
namespace Energy {

RegisteredPKFactory<AdvectionDiffusion> AdvectionDiffusion::reg_("advection-diffusion energy");

void AdvectionDiffusion::setup(const Teuchos::Ptr<State>& S) {

  // require fields
  Teuchos::RCP<CompositeVectorFactory> factory;

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
  factory->SetComponent("cell", AmanziMesh::CELL, 1);

  // boundary conditions
  Teuchos::ParameterList bc_plist = energy_plist_.sublist("boundary conditions", true);
  EnergyBCFactory bc_factory(S->GetMesh(), bc_plist);
  bc_temperature_ = bc_factory.CreateTemperature();
  bc_flux_ = bc_factory.CreateEnthalpyFlux();

  // operator for advection terms
  Operators::AdvectionFactory advection_factory;
  Teuchos::ParameterList advect_plist = energy_plist_.sublist("Advection");
  advection_ = advection_factory.create(advect_plist, S->GetMesh());
  advection_->set_num_dofs(1);

  // operator for the diffusion terms
  Teuchos::ParameterList mfd_plist = energy_plist_.sublist("Diffusion");
  matrix_ = Teuchos::rcp(new Operators::MatrixMFD(mfd_plist, S->GetMesh()));
  matrix_->SetSymmetryProperty(true);
  matrix_->SymbolicAssembleGlobalMatrices();
  matrix_->CreateMFDmassMatrices(Teuchos::null);

  // preconditioner
  // NOTE: may want to allow these to be the same/different?
  Teuchos::ParameterList mfd_pc_plist = energy_plist_.sublist("Diffusion PC");
  Teuchos::RCP<Operators::MatrixMFD> precon =
    Teuchos::rcp(new Operators::MatrixMFD(mfd_pc_plist, S->GetMesh()));
  set_preconditioner(precon);
};


// -- Initialize owned (dependent) variables.
void AdvectionDiffusion::initialize(const Teuchos::Ptr<State>& S) {
  PKPhysicalBDFBase::initialize(S);

  // initialize boundary conditions
  int nfaces = S->GetMesh()->num_entities(AmanziMesh::FACE, AmanziMesh::USED);
  bc_markers_.resize(nfaces, Operators::MATRIX_BC_NULL);
  bc_values_.resize(nfaces, 0.0);

  double time = S->time();
  bc_temperature_->Compute(time);
  bc_flux_->Compute(time);
  UpdateBoundaryConditions_();
};


// Evaluate BCs
void AdvectionDiffusion::UpdateBoundaryConditions_() {
  for (int n=0; n!=bc_markers_.size(); ++n) {
    bc_markers_[n] = Operators::MATRIX_BC_NULL;
    bc_values_[n] = 0.0;
  }

  Functions::BoundaryFunction::Iterator bc;
  for (bc=bc_temperature_->begin(); bc!=bc_temperature_->end(); ++bc) {
    int f = bc->first;
    bc_markers_[f] = Operators::MATRIX_BC_DIRICHLET;
    bc_values_[f] = bc->second;
  }

  for (bc=bc_flux_->begin(); bc!=bc_flux_->end(); ++bc) {
    int f = bc->first;
    bc_markers_[f] = Operators::MATRIX_BC_FLUX;
    bc_values_[f] = bc->second;
  }
};


// Add a boundary marker to owned faces.
void AdvectionDiffusion::ApplyBoundaryConditions_(const Teuchos::RCP<CompositeVector>& temperature) {
  int nfaces = temperature->size("face",false);
  for (int f=0; f!=nfaces; ++f) {
    if (bc_markers_[f] == Operators::MATRIX_BC_DIRICHLET) {
      (*temperature)("face", 0, f) = bc_values_[f];
    }
  }
};

} // namespace Energy
} // namespace Amanzi
