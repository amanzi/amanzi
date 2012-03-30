/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/*
This is the flow component of the Amanzi code. 
License: BSD
Authors: Neil Carlson (version 1) 
         Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
         Ethan Coon (ATS version) (ecoon@lanl.gov)
*/

#include "richards.hh"

namespace Amanzi {
namespace Flow {

// constructor
Richards::Richards(Teuchos::ParameterList& flow_plist, const Teuchos::RCP<State>& S,
           const Teuchos::RCP<TreeVector>& solution) :
    Flow(flow_plist,S,solution) {

  // just the extras...
  // data layouts for fields
  std::vector<AmanziMesh::Entity_kind> locations2(2);
  std::vector<std::string> names2(2);
  locations2[0] = AmanziMesh::CELL;
  locations2[1] = AmanziMesh::FACE;
  names2[0] = "cell";
  names2[1] = "face";

  // -- secondary variables
  S->RequireField("saturation_liquid", "flow", AmanziMesh::FACE, 1, true);
  S->RequireField("saturation_gas", "flow", AmanziMesh::FACE, 1, true);
  S->RequireField("density_gas", "flow", AmanziMesh::CELL, 1, true);
  S->RequireField("mol_frac_gas", "flow", AmanziMesh::CELL, 1, true);
  S->RequireField("molar_density_gas", "flow", AmanziMesh::CELL, 1, true);
  S->RequireField("molar_density_liquid", "flow", AmanziMesh::CELL, 1, true);
  S->RequireScalar("atmospheric_pressure");
  S->RequireField("relative_permeability", "flow", names2, locations2, 1, true);
};

// -- Initialize owned (dependent) variables.
void Richards::initialize(const Teuchos::RCP<State>& S) {

  // initial timestep size
  dt_ = flow_plist_.get<double>("Initial time step", 1.);

  // initialize pressure?
  if (flow_plist_.isParameter("Constant pressure")) {
    double p = flow_plist_.get<double>("Constant pressure");
    S->GetFieldData("pressure", "flow")->PutScalar(p);
    S->GetRecord("pressure", "flow")->set_initialized();
  }

  // update secondary variables for IC
  UpdateSecondaryVariables_(S);
  UpdatePermeability_(S);

  // declare secondary variables initialized
  S->GetRecord("saturation_liquid","flow")->set_intialized();
  S->GetRecord("saturation_gas","flow")->set_intialized();
  S->GetRecord("density_liquid","flow")->set_intialized();
  S->GetRecord("density_gas","flow")->set_intialized();
  S->GetRecord("molar_density_liquid","flow")->set_intialized();
  S->GetRecord("molar_density_gas","flow")->set_intialized();
  S->GetRecord("mol_frac_gas","flow")->set_intialized();
  S->GetRecord("relative_permeability","flow")->set_intialized();

  // initialize boundary conditions
  int nfaces = S->mesh()->num_entities(AmanziMesh::FACE, AmanziMesh::USED);
  bc_markers_.resize(nfaces, Operators::MFD_BC_NULL);
  bc_values_.resize(nfaces, 0.0);

  double time = S->time();
  bc_pressure_->Compute(time);
  bc_head_->Compute(time);
  bc_flux_->Compute(time);
  UpdateBoundaryConditions_();

  // initialize the timesteppper
  state_to_solution(S, solution_);
  if (!energy_plist_.get<bool>("Strongly Coupled PK", false)) {
    // model evaluator params
    // -- tolerances
    atol_ = flow_plist_.get<double>("Absolute error tolerance",1e-5);
    rtol_ = flow_plist_.get<double>("Relative error tolerance",1e-5);

    // -- instantiate time stepper
    Teuchos::RCP<Teuchos::ParameterList> bdf2_plist_p =
      Teuchos::rcp(new Teuchos::ParameterList(flow_plist_.sublist("Time integrator")));
    time_stepper_ = Teuchos::rcp(new BDF2TimeIntegrator(this, bdf2_plist_p, solution_));

    // -- initialize time derivative
    Teuchos::RCP<TreeVector> solution_dot = Teuchos::rcp( new TreeVector(*solution_));
    solution_dot->PutScalar(0.0);

    // -- set initial state
    time_stepper_->set_initial_state(S->time(), solution_, solution_dot);
  }
};

// -- Advance from state S0 to state S1 at time S0.time + dt.
bool Richards::advance(double dt) {
  state_to_solution(S_next_, solution_);

  // take the bdf timestep
  double h = dt;
  time_stepper_->time_step(h, solution_);
  time_stepper_->commit_solution(h, solution_);
  return false;
};
  
} // namespace
} // namespace

  
