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

AdvectionDiffusion::AdvectionDiffusion(Teuchos::ParameterList& plist,
        const Teuchos::RCP<State>& S,
        const Teuchos::RCP<TreeVector>& solution) :
    energy_plist_(plist) {

  solution_ = solution;

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
  factory->SetMesh(S->Mesh());
  factory->SetGhosted(true);
  factory->SetComponents(names2, locations2, num_dofs2);

  // -- thermal conductivity -- just cells
  factory = S->RequireField("thermal_conductivity", "energy");
  factory->SetMesh(S->Mesh());
  factory->SetGhosted(true);
  factory->SetComponent("cell", AmanziMesh::CELL, 1);

  // -- independent variables (not owned by this pk)
  factory = S->RequireField("darcy_flux");
  factory->SetMesh(S->Mesh());
  factory->SetGhosted(true);
  factory->SetComponent("face", AmanziMesh::FACE, 1);

  factory = S->RequireField("cell_volume");
  factory->SetMesh(S->Mesh());
  factory->SetGhosted(true);
  factory->SetComponent("cell", AmanziMesh::CELL, 1);

  // boundary conditions
  Teuchos::ParameterList bc_plist = energy_plist_.sublist("boundary conditions", true);
  EnergyBCFactory bc_factory(S->Mesh(), bc_plist);
  bc_temperature_ = bc_factory.CreateTemperature();
  bc_flux_ = bc_factory.CreateEnthalpyFlux();

  // operator for advection terms
  Operators::AdvectionFactory advection_factory;
  Teuchos::ParameterList advect_plist = energy_plist_.sublist("Advection");
  advection_ = advection_factory.create(advect_plist, S->Mesh());
  advection_->set_num_dofs(1);

  // operator for the diffusion terms
  Teuchos::ParameterList mfd_plist = energy_plist_.sublist("Diffusion");
  matrix_ = Teuchos::rcp(new Operators::MatrixMFD(mfd_plist, S->Mesh()));
  matrix_->SetSymmetryProperty(true);
  matrix_->SymbolicAssembleGlobalMatrices();

  // preconditioner
  // NOTE: may want to allow these to be the same/different?
  Teuchos::ParameterList mfd_pc_plist = energy_plist_.sublist("Diffusion PC");
  preconditioner_ = Teuchos::rcp(new Operators::MatrixMFD(mfd_pc_plist, S->Mesh()));
  preconditioner_->SetSymmetryProperty(true);
  preconditioner_->SymbolicAssembleGlobalMatrices();
  preconditioner_->InitPreconditioner(mfd_pc_plist);
};


// -- Initialize owned (dependent) variables.
void AdvectionDiffusion::initialize(const Teuchos::RCP<State>& S) {

  // initial timestep size
  dt_ = energy_plist_.get<double>("Initial time step", 1.);

  // initialize the operators
  matrix_->CreateMFDmassMatrices(Teuchos::null);
  preconditioner_->CreateMFDmassMatrices(Teuchos::null);

  // initialize boundary conditions
  int nfaces = S->Mesh()->num_entities(AmanziMesh::FACE, AmanziMesh::USED);
  bc_markers_.resize(nfaces, Operators::MFD_BC_NULL);
  bc_values_.resize(nfaces, 0.0);

  double time = S->time();
  bc_temperature_->Compute(time);
  bc_flux_->Compute(time);
  UpdateBoundaryConditions_();

  // update face temperature IC as a hint?
  Teuchos::RCP<CompositeVector> temp = S->GetFieldData("temperature", "energy");
  DeriveFaceValuesFromCellValues_(S, temp);

  // initialize time integrater and nonlinear solver
  solution_->set_data(temp);
  atol_ = energy_plist_.get<double>("Absolute error tolerance",1e0);
  rtol_ = energy_plist_.get<double>("Relative error tolerance",1e0);

  if (!energy_plist_.get<bool>("Strongly Coupled PK", false)) {
    // -- instantiate time stepper
    Teuchos::RCP<Teuchos::ParameterList> bdf2_plist_p =
      Teuchos::rcp(new Teuchos::ParameterList(energy_plist_.sublist("Time integrator")));
    time_stepper_ = Teuchos::rcp(new BDF2TimeIntegrator(this, bdf2_plist_p, solution_));

    // -- initialize time derivative
    Teuchos::RCP<TreeVector> solution_dot = Teuchos::rcp( new TreeVector(*solution_));
    solution_dot->PutScalar(0.0);

    // -- set initial state
    time_stepper_->set_initial_state(S->time(), solution_, solution_dot);
  }
};


// -- transfer operators -- ONLY COPIES POINTERS
void AdvectionDiffusion::state_to_solution(const Teuchos::RCP<State>& S,
        const Teuchos::RCP<TreeVector>& solution) {
  //Teuchos::RCP<CompositeVector> temp = S->GetFieldData("temperature", "energy");
  solution->set_data(S->GetFieldData("temperature", "energy"));
};

void AdvectionDiffusion::solution_to_state(const Teuchos::RCP<TreeVector>& solution,
        const Teuchos::RCP<State>& S) {
  //Teuchos::RCP<CompositeVector> temp = solution->data();
  S->SetData("temperature", "energy", solution->data());
};

  // -- Choose a time step compatible with physics.
double AdvectionDiffusion::get_dt() {
  return dt_;
};

// -- Advance from state S0 to state S1 at time S0.time + dt.
bool AdvectionDiffusion::advance(double dt) {
  state_to_solution(S_next_, solution_);

  // take the bdf timestep
  double h = dt;
  time_stepper_->time_step(h, solution_);

  // commit the step as successful
  time_stepper_->commit_solution(h, solution_);
  solution_to_state(solution_, S_next_);
  commit_state(h, S_next_);

  return false;
};


// -----------------------------------------------------------------------------
// Interpolate pressure ICs on cells to ICs for lambda (faces).
// -----------------------------------------------------------------------------
void AdvectionDiffusion::DeriveFaceValuesFromCellValues_(const Teuchos::RCP<State>& S,
        const Teuchos::RCP<CompositeVector>& temp) {
  AmanziMesh::Entity_ID_List cells;

  int f_owned = temp->size("face");
  for (int f=0; f!=f_owned; ++f) {
    cells.clear();
    S->Mesh()->face_get_cells(f, AmanziMesh::USED, &cells);
    int ncells = cells.size();

    double face_value = 0.0;
    for (int n=0; n!=ncells; ++n) {
      face_value += (*temp)("cell",cells[n]);
    }
    (*temp)("face",f) = face_value / ncells;
  }
};


// Evaluate BCs
void AdvectionDiffusion::UpdateBoundaryConditions_() {
  for (int n=0; n!=bc_markers_.size(); ++n) {
    bc_markers_[n] = Operators::MFD_BC_NULL;
    bc_values_[n] = 0.0;
  }

  Functions::BoundaryFunction::Iterator bc;
  for (bc=bc_temperature_->begin(); bc!=bc_temperature_->end(); ++bc) {
    int f = bc->first;
    bc_markers_[f] = Operators::MFD_BC_DIRICHLET;
    bc_values_[f] = bc->second;
  }

  for (bc=bc_flux_->begin(); bc!=bc_flux_->end(); ++bc) {
    int f = bc->first;
    bc_markers_[f] = Operators::MFD_BC_FLUX;
    bc_values_[f] = bc->second;
  }
};


// Add a boundary marker to owned faces.
void AdvectionDiffusion::ApplyBoundaryConditions_(const Teuchos::RCP<CompositeVector>& temperature) {
  int nfaces = temperature->size("face",false);
  for (int f=0; f!=nfaces; ++f) {
    if (bc_markers_[f] == Operators::MFD_BC_DIRICHLET) {
      (*temperature)("face", 0, f) = bc_values_[f];
    }
  }
};

} // namespace Energy
} // namespace Amanzi
