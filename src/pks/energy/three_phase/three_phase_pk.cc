/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon
------------------------------------------------------------------------- */

#include "bdf1_time_integrator.hh"
#include "advection_factory.hh"
#include "energy_bc_factory.hh"
#include "thermal_conductivity_threephase_factory.hh"
#include "iem_factory.hh"

#include "three_phase.hh"

namespace Amanzi {
namespace Energy {

RegisteredPKFactory<ThreePhase> ThreePhase::reg_("three-phase energy");

// -------------------------------------------------------------
// Constructor
// -------------------------------------------------------------
ThreePhase::ThreePhase(Teuchos::ParameterList& plist,
                       const Teuchos::RCP<State>& S,
                       const Teuchos::RCP<TreeVector>& solution) :
    energy_plist_(plist) {

  solution_ = solution;

  // require fields
  // primary variable: temperature on both cells and faces, ghosted, with 1 dof
  std::vector<AmanziMesh::Entity_kind> locations2(2);
  std::vector<std::string> names2(2);
  std::vector<int> num_dofs2(2,1); // = [1, 1]
  locations2[0] = AmanziMesh::CELL;
  locations2[1] = AmanziMesh::FACE;
  names2[0] = "cell";
  names2[1] = "face";

  S->RequireField("temperature", "energy")->SetMesh(S->Mesh())->SetGhosted()
                    ->SetComponents(names2, locations2, num_dofs2);

  // secondary variables
  //   Nearly all of these will be defined on a single dof on cells.  To make
  //   this process easier, we set up two factories, one owned, one not, and
  //   use those to initialize the data factories.

  CompositeVectorFactory one_cell_owned_factory;
  one_cell_owned_factory.SetMesh(S->Mesh());
  one_cell_owned_factory.SetGhosted();
  one_cell_owned_factory.SetComponent("cell", AmanziMesh::CELL, 1);
                          // The use of SetComponent() implies this is final.

  CompositeVectorFactory one_cell_factory;
  one_cell_owned_factory.SetMesh(S->Mesh());
  one_cell_owned_factory.SetGhosted();
  one_cell_owned_factory.AddComponent("cell", AmanziMesh::CELL, 1);
                          // The use of AddComponent() implies the actual
                          // vector may include other components.

  // -- gas
  S->RequireField("internal_energy_gas", "energy")->Update(one_cell_owned_factory);

  // -- liquid
  S->RequireField("internal_energy_liquid", "energy")->Update(one_cell_owned_factory);
  S->RequireField("enthalpy_liquid", "energy")->Update(one_cell_owned_factory);

  // -- ice
  S->RequireField("internal_energy_ice", "energy")->Update(one_cell_owned_factory);

  // -- rock assumed constant for now?
  S->RequireScalar("density_rock");
  S->RequireField("internal_energy_rock", "energy")->Update(one_cell_owned_factory);

  // -- thermal conductivity
  S->RequireField("thermal_conductivity", "energy")->Update(one_cell_owned_factory);

  // independent variables
  //   These are not owned by this pk, but we want to make sure they exist,
  //   on our mesh, on cells.
  S->RequireField("porosity")->Update(one_cell_factory);
  S->RequireField("density_gas")->Update(one_cell_factory);
  S->RequireField("density_liquid")->Update(one_cell_factory);
  S->RequireField("density_ice")->Update(one_cell_factory);
  S->RequireField("molar_density_gas")->Update(one_cell_factory);
  S->RequireField("molar_density_liquid")->Update(one_cell_factory);
  S->RequireField("molar_density_ice")->Update(one_cell_factory);
  S->RequireField("mol_frac_gas")->Update(one_cell_factory);
  S->RequireField("saturation_gas")->Update(one_cell_factory);
  S->RequireField("saturation_liquid")->Update(one_cell_factory);
  S->RequireField("saturation_ice")->Update(one_cell_factory);

  S->RequireField("darcy_flux")->SetMesh(S->Mesh())->SetGhosted()
                                ->AddComponent("face", AmanziMesh::FACE, 1);

  S->RequireField("pressure")->SetMesh(S->Mesh())->SetGhosted()
                                ->AddComponent("cell", AmanziMesh::CELL, 1);

  // abs conductivity tensor
  int c_owned = S->Mesh()->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);
  Ke_.resize(c_owned);
  for (int c=0; c!=c_owned; ++c) {
    Ke_[c].init(S->Mesh()->space_dimension(),1);
    Ke_[c](1,1) = 1.0;
  }

  // constitutive relations
  // -- thermal conductivity
  Teuchos::ParameterList tcm_plist = energy_plist_.sublist("Thermal Conductivity Model");
  EnergyRelations::ThermalConductivityThreePhaseFactory tc_factory;
  thermal_conductivity_model_ = tc_factory.createThermalConductivityModel(tcm_plist);

  // -- internal energy models
  // -- water vapor requires a specialized model
  Teuchos::ParameterList ieg_plist = energy_plist_.sublist("Internal Energy Gas Model");
  iem_gas_ =
    Teuchos::rcp(new EnergyRelations::InternalEnergyWaterVapor(ieg_plist));

  // -- other IEM come from factory
  // -- for liquid
  EnergyRelations::IEMFactory iem_factory;
  Teuchos::ParameterList iel_plist = energy_plist_.sublist("Internal Energy Liquid Model");
  iem_liquid_ = iem_factory.createIEM(iel_plist);

  // -- for ice
  Teuchos::ParameterList iei_plist = energy_plist_.sublist("Internal Energy Ice Model");
  iem_ice_ = iem_factory.createIEM(iei_plist);

  // -- for rock
  Teuchos::ParameterList ier_plist = energy_plist_.sublist("Internal Energy Rock Model");
  iem_rock_ = iem_factory.createIEM(ier_plist);

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
  bool symmetric = true;
  Teuchos::ParameterList mfd_plist = energy_plist_.sublist("Diffusion");
  matrix_ = Teuchos::rcp(new Operators::MatrixMFD(mfd_plist, S->Mesh()));
  matrix_->SetSymmetryProperty(symmetric);
  matrix_->SymbolicAssembleGlobalMatrices();
  matrix_->CreateMFDmassMatrices(Ke_);

  // preconditioner
  // NOTE: may want to allow these to be the same/different?
  Teuchos::ParameterList mfd_pc_plist = energy_plist_.sublist("Diffusion PC");
  preconditioner_ = Teuchos::rcp(new Operators::MatrixMFD(mfd_pc_plist, S->Mesh()));
  preconditioner_->SetSymmetryProperty(symmetric);
  preconditioner_->SymbolicAssembleGlobalMatrices();
  preconditioner_->CreateMFDmassMatrices(Ke_);

  Teuchos::ParameterList mfd_pc_ml_plist = mfd_pc_plist.sublist("ML Parameters");
  preconditioner_->InitMLPreconditioner(mfd_pc_ml_plist);
};


// -------------------------------------------------------------
// Initialize PK
// -------------------------------------------------------------
void ThreePhase::initialize(const Teuchos::RCP<State>& S) {

  // initial timestep size
  dt_ = energy_plist_.get<double>("Initial time step", 1.);

  // initialize boundary conditions
  int nfaces = S->Mesh()->num_entities(AmanziMesh::FACE, AmanziMesh::USED);
  bc_markers_.resize(nfaces, Operators::MFD_BC_NULL);
  bc_values_.resize(nfaces, 0.0);

  // update face temperatures as a hint?
  Teuchos::RCP<CompositeVector> temp = S->GetFieldData("temperature", "energy");
  DeriveFaceValuesFromCellValues_(S, temp);

  // declare secondary variables initialized, as they get done in a call to
  // commit_state
  S->GetField("thermal_conductivity","energy")->set_initialized();
  S->GetField("internal_energy_gas","energy")->set_initialized();
  S->GetField("internal_energy_liquid","energy")->set_initialized();
  S->GetField("internal_energy_ice","energy")->set_initialized();
  S->GetField("internal_energy_rock","energy")->set_initialized();
  S->GetField("enthalpy_liquid","energy")->set_initialized();

  // initialize the timesteppper
  solution_->set_data(temp);
  atol_ = energy_plist_.get<double>("Absolute error tolerance",1e0);
  rtol_ = energy_plist_.get<double>("Relative error tolerance",1e0);

  if (!energy_plist_.get<bool>("Strongly Coupled PK", false)) {
    // -- instantiate time stepper
    Teuchos::RCP<Teuchos::ParameterList> bdf1_plist_p =
      Teuchos::rcp(new Teuchos::ParameterList(energy_plist_.sublist("Time integrator")));
    time_stepper_ = Teuchos::rcp(new BDF1TimeIntegrator(this, bdf1_plist_p, solution_));
    time_step_reduction_factor_ = bdf1_plist_p->get<double>("time step reduction factor");

    // -- initialize time derivative
    Teuchos::RCP<TreeVector> solution_dot = Teuchos::rcp( new TreeVector(*solution_));
    solution_dot->PutScalar(0.0);

    // -- set initial state
    time_stepper_->set_initial_state(S->time(), solution_, solution_dot);
  }
};


// -----------------------------------------------------------------------------
// Update any secondary (dependent) variables given a solution.
//
//   After a timestep is evaluated (or at ICs), there is no way of knowing if
//   secondary variables have been updated to be consistent with the new
//   solution.
// -----------------------------------------------------------------------------
void ThreePhase::commit_state(double dt, const Teuchos::RCP<State>& S) {
  // update secondary variables for IC i/o
  UpdateSecondaryVariables_(S);
  UpdateThermalConductivity_(S);
  UpdateEnthalpyLiquid_(S);
};


// -----------------------------------------------------------------------------
// Transfer operators -- ONLY COPIES POINTERS
// -----------------------------------------------------------------------------
void ThreePhase::state_to_solution(const Teuchos::RCP<State>& S,
        const Teuchos::RCP<TreeVector>& solution) {
  solution->set_data(S->GetFieldData("temperature", "energy"));
};


// -----------------------------------------------------------------------------
// Transfer operators -- ONLY COPIES POINTERS
// -----------------------------------------------------------------------------
void ThreePhase::solution_to_state(const Teuchos::RCP<TreeVector>& solution,
        const Teuchos::RCP<State>& S) {
  S->SetData("temperature", "energy", solution->data());
};


// -----------------------------------------------------------------------------
// Choose a time step compatible with physics.
// -----------------------------------------------------------------------------
double ThreePhase::get_dt() {
  return dt_;
};


// -----------------------------------------------------------------------------
// Advance from state S to state S_next at time S.time + dt.
// -----------------------------------------------------------------------------
bool ThreePhase::advance(double dt) {
  state_to_solution(S_next_, solution_);

  // take a bdf timestep
  double h = dt;
  try {
    dt_ = time_stepper_->time_step(h, solution_);
  } catch (Exceptions::Amanzi_exception &error) {
    if (error.what() == std::string("BDF time step failed")) {
      // try cutting the timestep
      dt_ = dt_*time_step_reduction_factor_;
      return true;
    } else {
      throw error;
    }
  }

  // commit the step as successful
  time_stepper_->commit_solution(h, solution_);
  solution_to_state(solution_, S_next_);
  commit_state(h, S_next_);

  return false;
};


// -----------------------------------------------------------------------------
// Interpolate pressure ICs on cells to ICs for lambda (faces).
// -----------------------------------------------------------------------------
void ThreePhase::DeriveFaceValuesFromCellValues_(const Teuchos::RCP<State>& S,
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


// -----------------------------------------------------------------------------
// Evaluate boundary conditions at the current time.
// -----------------------------------------------------------------------------
void ThreePhase::UpdateBoundaryConditions_() {
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


// -----------------------------------------------------------------------------
// Add a boundary marker to owned faces.
// -----------------------------------------------------------------------------
void ThreePhase::ApplyBoundaryConditions_(const Teuchos::RCP<CompositeVector>& temperature) {
  int nfaces = temperature->size("face");
  for (int f=0; f!=nfaces; ++f) {
    if (bc_markers_[f] == Operators::MFD_BC_DIRICHLET) {
      (*temperature)("face",f) = bc_values_[f];
    }
  }
};

} // namespace Energy
} // namespace Amanzi
