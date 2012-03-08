/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon
------------------------------------------------------------------------- */

#include "bdf2_time_integrator.hh"
#include "advection_factory.hh"

#include "air_water_rock.hh"

namespace Amanzi {
namespace Energy {


AirWaterRock::AirWaterRock(Teuchos::ParameterList& plist,
        Teuchos::RCP<State>& S, Teuchos::RCP<TreeVector>& soln) :
    energy_plist_(plist) {

  solution_ = soln;

  // require fields
  // primary variable: temperature on both cells and faces, ghosted, with 1 dof
  std::vector<AmanziMesh::Entity_kind> locations(2);
  std::vector<std::string> names(2);
  std::vector< std::vector<std::string> > subfield_names(2);
  locations[0] = AmanziMesh::CELL;
  locations[1] = AmanziMesh::FACE;
  names[0] = "cell";
  names[1] = "face";
  subfield_names[0].resize(1); subfield_names[0][0] = "temperature";
  subfield_names[1].resize(1); subfield_names[1][0] = "temperature_lambda";

  S->RequireField("temperature", "energy", names, locations, 1, true);
  S->GetRecord("temperature","energy")->set_io_vis(true);
  Teuchos::RCP<CompositeVector> temp = S->GetFieldData("temperature", "energy");
  temp->set_subfield_names(subfield_names);
  solution_->set_data(temp);

  // secondary variables
  // -- gas
  S->RequireField("internal_energy_gas", "energy", AmanziMesh::CELL, 1, true);

  // -- liquid
  S->RequireField("internal_energy_liquid", "energy", AmanziMesh::CELL, 1, true);
  S->RequireField("specific_enthalpy_liquid", "energy", AmanziMesh::CELL, 1, true);

  // -- rock assumed constant for now?
  S->RequireScalar("density_rock");
  S->RequireField("internal_energy_rock", "energy", AmanziMesh::CELL, 1, true);


  // -- parameters
  // S->RequireField("thermal_conductivity", "energy", names, locations, 1, true); // if upwinded
  S->RequireField("thermal_conductivity", "energy", AmanziMesh::CELL, 1, true); // if not

  // independent variables (not owned by this pk)
  S->RequireField("porosity", AmanziMesh::CELL, 1, true);
  S->RequireField("density_gas", AmanziMesh::CELL, 1, true);
  S->RequireField("density_liquid", AmanziMesh::CELL, 1, true);
  S->RequireField("saturation_liquid", AmanziMesh::CELL, 1, true);
  S->RequireField("darcy_flux", AmanziMesh::FACE, 1, true);


  // operator for advection terms
  Operators::AdvectionFactory advection_factory;
  Teuchos::ParameterList advect_plist = energy_plist_.sublist("Advection");
  advection_ = advection_factory.create(advect_plist, S->mesh());
  advection_->set_num_dofs(1);

  // check if we need to make a time integrator
  if (!energy_plist_.get<bool>("Strongly Coupled PK", false)) {
    Teuchos::RCP<Teuchos::ParameterList> bdf2_plist_p =
      Teuchos::rcp(new Teuchos::ParameterList(energy_plist_.sublist("Time integrator")));
    time_stepper_ = Teuchos::rcp(new BDF2TimeIntegrator(this, bdf2_plist_p, solution_));
  }
};

// -- Initialize owned (dependent) variables.
void AirWaterRock::initialize(const Teuchos::RCP<State>& S) {
  // constant initial temperature
  if (energy_plist_.isParameter("Constant temperature")) {
    double T = energy_plist_.get<double>("Constant temperature");
    S->GetFieldData("temperature", "energy")->PutScalar(T);
    S->GetRecord("temperature", "energy")->set_initialized();
  }
  dt_ = energy_plist_.get<double>("time step size",1.0);

  state_to_solution(S, solution_);

  if (!energy_plist_.get<bool>("Strongly Coupled PK", false)) {
    // model evaluator params
    // -- tolerances
    atol_ = energy_plist_.get<double>("Absolute error tolerance",1e-5);
    rtol_ = energy_plist_.get<double>("Relative error tolerance",1e-5);

    // -- initialize time derivative
    Teuchos::RCP<TreeVector> solution_dot = Teuchos::rcp( new TreeVector(*solution_));
    solution_dot->PutScalar(0.0);

    // -- set initial state
    time_stepper_->set_initial_state(S->time(), solution_, solution_dot);
  }
};

// -- transfer operators -- ONLY COPIES POINTERS
void AirWaterRock::state_to_solution(const Teuchos::RCP<State>& S,
        const Teuchos::RCP<TreeVector>& soln) {
  //Teuchos::RCP<CompositeVector> temp = S->GetFieldData("temperature", "energy");
  soln->set_data(S->GetFieldData("temperature", "energy"));
};

void AirWaterRock::solution_to_state(const Teuchos::RCP<TreeVector>& soln,
        const Teuchos::RCP<State>& S) {
  //Teuchos::RCP<CompositeVector> temp = soln->data();
  S->SetData("temperature", "energy", soln->data());
};

  // -- Choose a time step compatible with physics.
double AirWaterRock::get_dt() {
  return dt_;
};

// -- Advance from state S0 to state S1 at time S0.time + dt.
bool AirWaterRock::advance(double dt) {
  // take the bdf timestep
  double h = dt;
  time_stepper_->time_step(h, solution_);
  time_stepper_->commit_solution(h, solution_);
  return false;
};

} // namespace Energy
} // namespace Amanzi
