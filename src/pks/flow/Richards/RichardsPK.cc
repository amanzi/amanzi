/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
   ATS

   License: see $ATS_DIR/COPYRIGHT
   Author: Ethan Coon

   Implementation for a basic Richards PK.

   Example usage:

   <ParameterList name="flow">
   <Parameter name="PK model" type="string" value="Richards"/>
   </ParameterList>

   ------------------------------------------------------------------------- */

#include "RichardsPK.hh"
#include "RichardsProblem.hh"

namespace Amanzi {


/* ******************************************************************
*  Constructor for the basic Richards Process Kernel
*********************************************************************/
RichardsPK::RichardsPK(Teuchos::ParameterList &flow_plist, Teuchos::RCP<State> &S,
                       Teuchos::RCP<TreeVector>& solution) :
    flow_plist_(flow_plist) {

  solution_ = solution;

  // Require fields in the state
  // -- pressure field
  S->require_field("pressure", Amanzi::AmanziMesh::CELL, "flow");
  S->get_field_record("pressure")->set_io_vis(true);
  Teuchos::RCP<Epetra_MultiVector> pressure_ptr = S->get_field("pressure", "flow");
  solution_->PushBack(pressure_ptr);

  // -- constraints -- pressure on the faces
  S->require_field("pressure_lambda", Amanzi::AmanziMesh::FACE, "flow");
  Teuchos::RCP<Epetra_MultiVector> pressure_lambda_ptr =
    S->get_field("pressure_labmda", "flow");
  solution_->PushBack(pressure_lambda_ptr);

  // -- liquid flux, scalar on faces
  S->require_field("darcy_flux", Amanzi::AmanziMesh::FACE, "flow");

  // -- darcy velocity, vector on cells (for diagnostics only)
  S->require_field("darcy_velocity", Amanzi::AmanziMesh::CELL, "flow", 3);
  S->get_field_record("darcy_velocity")->set_io_vis(true);
  std::vector<std::string> subfieldnames;
  subfieldnames.push_back("Liquid X-Velocity");
  subfieldnames.push_back("Liquid Y-Velocity");
  subfieldnames.push_back("Liquid Z-Velocity");
  S->set_subfield_names("darcy_velocity", subfieldnames);

  // -- saturation
  S->require_field("saturation", Amanzi::AmanziMesh::CELL, "flow");
  S->get_field_record("saturation")->set_io_vis(true);

  // -- parameters
  S->require_field("permeability", Amanzi::AmanziMesh::CELL);
  S->require_field("relative_permeability", Amanzi::AmanziMesh::CELL);
  S->require_field("porosity", Amanzi::AmanziMesh::CELL);

  // create the problem
  Teuchos::ParameterList richards_plist = plist.sublist("Richards Problem");
  problem_ = Teuchos::rcp(new RichardsProblem(S->get_mesh_maps(), richards_plist));

  // check if we need to make a time integrator
  if (!flow_plist_.get<bool>("Strongly Coupled PK", false)) {
    time_stepper_ = Teuchos::rcp(new ImplicitTIBDF2(*this,solution_));
    Teuchos::RCP<Teuchos::ParameterList> bdf2_list_p(new Teuchos::ParameterList(flow_plist_.sublist("Time integrator")));
    time_stepper_->setParameterList(bdf2_list_p);
  }
};

/* ******************************************************************
*  Initialization of state.
*********************************************************************/
void RichardsPK::initialize(Teuchos::RCP<State> &S) {
  Teuchos::ParameterList richards_plist = plist.sublist("Richards Problem");

  // set the timestep sizes
  h0_ = richards_plist.get<double>("Initial time step");
  h_ = h0_;
  hnext_ = h0_;

  // Set the problem independent variables from state
  problem_->SetFluidDensity(*S->get_density());
  problem_->SetFluidViscosity(*S->get_viscosity());
  problem_->SetGravity(*S->get_gravity());
  problem_->Initialize();

  // Initialize the problem, and calculate initial values for dependent
  // variables.  In the process, set the fields to initialized.
  // -- IC for pressure
  Teuchos::RCP<Epetra_MultiVector> pressure = S->get_field("pressure", "flow");
  Teuchos::RCP<Epetra_MultiVector> pressure_lambda = S->get_field("pressure_lambda", "flow");
  problem_->InitializePressureCells(pressure);
  problem_->InitializePressureFaces(pressure_lambda);
  S->get_field_record("pressure")->set_initialized();
  S->get_field_record("pressure_lambda")->set_initialized();

  // -- IC for saturations, if not set in state
  Teuchos::RCP<Epetra_MultiVector> saturation = S->get_field("saturation", "flow");
  problem_->InitializeSaturation(saturation);
  S->get_field_record("saturation")->set_initialized();

  // -- calculate rel perm
  Teuchos::RCP<Epetra_MultiVector> rel_perm = S->get_field("permeability", "flow");
  problem_->DeriveRelPerm(*saturation, rel_perm);
  S->get_field_record("relative_permeability")->set_initialized();

  // -- calculate darcy flux
  ASSERT(S->get_field_record("permeability")->initialized());
  Teuchos::RCP<Epetra_MultiVector> perm = S->get_field("permeability");
  Teuchos::RCP<Epetra_MultiVector> darcy_flux = S->get_field("darcy_flux", "flow");
  // TODO: scoped pointer for l1_error
  double l1_error;
  problem_->DeriveDarcyFlux(*pressure, *rel_perm, *perm, darcy_flux, &l1_error);
  S->get_field_record("darcy_flux")->set_initialized();

  // -- darcy velocity
  Teuchos::RCP<Epetra_MultiVector> darcy_velocity = S->get_field("darcy_velocity", "flow");
  // TODO: scoped pointer for l1_error
  problem_->DeriveDarcyVelocity(*darcy_flux, darcy_velocity);
  S->get_field_record("darcy_velocity")->set_initialized();


  ///// STOP HERE ////////



  
  // set up preconditioner
  int errc;
  double t0 = S->get_time();
  problem->Compute_udot(t0, *solution, udot);
  time_stepper->set_initial_state(t0, *solution, udot);

  update_precon(t0, *solution, h, errc);

  // if steady state solve
  if (steady_state) {
    // NOTE: this does NOT get (nor need) state, to gaurantee it cannot alter State
    advance_to_steady_state();

    // update state
    S->set_field("pressure", "flow", *pressure_cells);
    S->set_field("darcy_flux", "flow", *richards_flux);
    S->set_time(ss_t1);
  }

  // IC for velocity
  Teuchos::RCP<Epetra_MultiVector> velocity = S->get_field("darcy_velocity", "flow");
  problem->DeriveDarcyVelocity(*solution, *velocity);
};

void RichardsPK::advance_to_steady_state()
{
  double t0 = ss_t0;
  double t1 = ss_t1;
  double h = h0;
  double hnext;

  // iterate
  int i = 0;
  double tlast = t0;

  do {
    time_stepper->bdf2_step(h,0.0,*solution,hnext);
    time_stepper->commit_solution(h,*solution);
    time_stepper->write_bdf2_stepping_statistics();

    h = hnext;
    i++;
    tlast=time_stepper->most_recent_time();
  } while (t1 >= tlast);

  // Derive the Richards fluxes on faces
  double l1_error;
  problem->DeriveDarcyFlux(*solution, *richards_flux, l1_error);
  std::cout << "L1 norm of the Richards flux discrepancy = " << l1_error << std::endl;
}

// NOTE: this should get split into three parts --
//    1. set problem data from state S0
//    2. advance transient (which in no way can alter either S0 or S1)
//    3. set state S1 with solution from problem
// much like advance_to_steady_state() works
bool RichardsPK::advance_transient(double dt, const Teuchos::RCP<State> &S0,
                                              Teuchos::RCP<State> &S1) {
  // potentially changed permeability and porosity
  problem->SetPermeability(*(*(S0->get_field("permeability")))(0));
  problem->SetPorosity(*(*(S0->get_field("porosity")))(0));

  // take the timestep
  time_stepper->bdf2_step(dt, 0.0, *solution, hnext);
  time_stepper->commit_solution(dt, *solution);
  time_stepper->write_bdf2_stepping_statistics();

  // update the state at the new time
  // update pressure
  S1->set_field("pressure", "flow", *pressure_cells);

  // update flux
  double l1_error;
  problem->DeriveDarcyFlux(*solution, *richards_flux, l1_error);
  std::cout << "L1 norm of the Richards flux discrepancy = " << l1_error << std::endl;
  S1->set_field("darcy_flux", "flow", *richards_flux);

  // update velocity DOES THIS NEED TO BE DONE NOW?  save for commit_state?
  Teuchos::RCP<Epetra_MultiVector> velocity = S1->get_field("darcy_velocity", "flow");
  problem->DeriveDarcyVelocity(*solution, *velocity);

  // update saturation
  Teuchos::RCP<Epetra_MultiVector> saturation = S1->get_field("water saturation", "flow");
  problem->DeriveSaturation(*pressure_cells, *(*saturation)(0));

  return false; // bdf2_step subcycles, so we have always succeeded
};
} // close namespace Amanzi
