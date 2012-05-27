/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon
------------------------------------------------------------------------- */

#include "bdf2_time_integrator.hh"
#include "advection_factory.hh"
#include "energy_bc_factory.hh"

#include "advection_diffusion.hh"

namespace Amanzi {
namespace Energy {

RegisteredPKFactory<AdvectionDiffusion> AdvectionDiffusion::reg_("advection-diffusion energy");

AdvectionDiffusion::AdvectionDiffusion(Teuchos::ParameterList& plist,
        const Teuchos::RCP<State>& S, const Teuchos::RCP<TreeVector>& solution) :
    energy_plist_(plist) {

  // require fields
  // primary variable: temperature on both cells and faces, ghosted, with 1 dof
  std::vector<AmanziMesh::Entity_kind> locations2(2);
  std::vector<std::string> names2(2);
  std::vector< std::vector<std::string> > subfield_names(2);
  locations2[0] = AmanziMesh::CELL;
  locations2[1] = AmanziMesh::FACE;
  names2[0] = "cell";
  names2[1] = "face";
  subfield_names[0].resize(1); subfield_names[0][0] = "temperature";
  subfield_names[1].resize(1); subfield_names[1][0] = "temperature_lambda";

  S->RequireField("temperature", "energy", names2, locations2, 1, true);
  S->GetRecord("temperature","energy")->set_io_vis(true);
  Teuchos::RCP<CompositeVector> temp = S->GetFieldData("temperature", "energy");
  temp->set_subfield_names(subfield_names);
  solution->set_data(temp);
  solution_ = solution;

  // -- parameters
  S->RequireField("thermal_conductivity", "energy", AmanziMesh::CELL, 1, true);

  // independent variables (not owned by this pk)
  S->RequireField("darcy_flux", AmanziMesh::FACE, 1, true);
  S->RequireField("cell_volume",AmanziMesh::CELL,1,true);

  // boundary conditions
  Teuchos::ParameterList bc_plist = energy_plist_.sublist("boundary conditions", true);
  EnergyBCFactory bc_factory(S->mesh(), bc_plist);
  bc_temperature_ = bc_factory.CreateTemperature();
  bc_flux_ = bc_factory.CreateEnthalpyFlux();

  // operator for advection terms
  Operators::AdvectionFactory advection_factory;
  Teuchos::ParameterList advect_plist = energy_plist_.sublist("Advection");
  advection_ = advection_factory.create(advect_plist, S->mesh());
  advection_->set_num_dofs(1);

  // operator for the diffusion terms
  Teuchos::ParameterList mfd_plist = energy_plist_.sublist("Diffusion");
  matrix_ = Teuchos::rcp(new Operators::MatrixMFD(mfd_plist, S->mesh()));
  matrix_->SetSymmetryProperty(true);
  matrix_->SymbolicAssembleGlobalMatrices();

  // preconditioner
  // NOTE: may want to allow these to be the same/different?
  Teuchos::ParameterList mfd_pc_plist = energy_plist_.sublist("Diffusion PC");
  preconditioner_ = Teuchos::rcp(new Operators::MatrixMFD(mfd_pc_plist, S->mesh()));
  preconditioner_->SetSymmetryProperty(true);
  preconditioner_->SymbolicAssembleGlobalMatrices();
  Teuchos::ParameterList mfd_pc_ml_plist = mfd_pc_plist.sublist("ML Parameters");
  preconditioner_->InitMLPreconditioner(mfd_pc_ml_plist);
};

// -- Initialize owned (dependent) variables.
void AdvectionDiffusion::initialize(const Teuchos::RCP<State>& S) {

  // initial timestep size
  dt_ = energy_plist_.get<double>("Initial time step", 1.);

  // constant initial temperature
  if (energy_plist_.isParameter("Constant temperature")) {
    double T = energy_plist_.get<double>("Constant temperature");
    S->GetFieldData("temperature", "energy")->PutScalar(T);
    S->GetRecord("temperature", "energy")->set_initialized();
  }

  // initialize thermal conductivity
  Teuchos::RCP<CompositeVector> thermal_conductivity =
    S->GetFieldData("thermal_conductivity", "energy");
  if (energy_plist_.isParameter("Constant thermal_conductivity")) {
    double K = energy_plist_.get<double>("Constant thermal_conductivity");
    thermal_conductivity->ViewComponent("cell")->PutScalar(K);
    S->GetRecord("thermal_conductivity", "energy")->set_initialized();
  }
  int size = thermal_conductivity->ViewComponent("cell",false)->MyLength();
  Ke_.resize(size);
  for (int c=0; c!=size; ++c) {
    Ke_[c].init(S->mesh()->space_dimension(),1);
  }

  // initialize boundary conditions
  int nfaces = S->mesh()->num_entities(AmanziMesh::FACE, AmanziMesh::USED);
  bc_markers_.resize(nfaces, Operators::MFD_BC_NULL);
  bc_values_.resize(nfaces, 0.0);

  double time = S->time();
  bc_temperature_->Compute(time);
  bc_flux_->Compute(time);
  UpdateBoundaryConditions_();

  state_to_solution(S, solution_);
  atol_ = energy_plist_.get<double>("Absolute error tolerance",1e-5);
  rtol_ = energy_plist_.get<double>("Relative error tolerance",1e-5);

  // initialize the timesteppper
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
  time_stepper_->commit_solution(h, solution_);
  return false;
};

void AdvectionDiffusion::UpdateBoundaryConditions_() {
  for (int n=0; n!=bc_markers_.size(); ++n) {
    bc_markers_[n] = Operators::MFD_BC_NULL;
    bc_values_[n] = 0.0;
  }

  BoundaryFunction::Iterator bc;
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

void AdvectionDiffusion::ApplyBoundaryConditions_(const Teuchos::RCP<CompositeVector>& temperature) {
  int nfaces = S_next_->mesh()->num_entities(AmanziMesh::FACE, AmanziMesh::OWNED);
  for (int f=0; f!=nfaces; ++f) {
    if (bc_markers_[f] == Operators::MFD_BC_DIRICHLET) {
      (*temperature)("face", 0, f) = bc_values_[f];
    }
  }
};

} // namespace Energy
} // namespace Amanzi
