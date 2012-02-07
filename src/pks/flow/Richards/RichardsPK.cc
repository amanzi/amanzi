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
  S->require_field("pressure_dot", Amanzi::AmanziMesh::CELL, "flow");

  // -- constraints -- pressure on the faces
  S->require_field("pressure_lambda", Amanzi::AmanziMesh::FACE, "flow");
  Teuchos::RCP<Epetra_MultiVector> pressure_lambda_ptr =
    S->get_field("pressure_labmda", "flow");
  solution_->PushBack(pressure_lambda_ptr);
  S->require_field("pressure_lambda_dot", Amanzi::AmanziMesh::FACE, "flow");

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
    time_stepper_ = Teuchos::rcp(new ImplicitTIBDF2(*this, solution_));
    Teuchos::RCP<Teuchos::ParameterList> bdf2_list_p(new Teuchos::ParameterList(flow_plist_.sublist("Time integrator")));
    time_stepper_->setParameterList(bdf2_list_p);
  }
};

/* ******************************************************************
*  Initialization of the problem, time stepper, and state.
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

  // set up timestepper
  Teuchos::RCP<TreeVector> solution_dot = Teuchos::rcp(new TreeVector(*solution_));
  state_to_solution(S, solution_, solution_dot);
  time_stepper_->set_initial_state(S->get_time(), solution_, solution_dot);

  // set up preconditioner
  int errc;
  // TODO: scoped pointer for errc
  update_precon(S->get_time(), solution_, h_, &errc);
};

  /* ******************************************************************
   *  Pointer copy of state to solution
   *********************************************************************/
  void RichardsPK::state_to_solution(Teuchos::RCP<State>& S,
                                     Teuchos::RCP<TreeVector>& soln) {
    (*soln)[0] = S->get_field("pressure", "flow");
    (*soln)[1] = S->get_field("pressure_lambda", "flow");
  };

  /* ******************************************************************
   *  Pointer copy of state to solution
   *********************************************************************/
  void RichardsPK::state_to_solution(Teuchos::RCP<State>& S,
                                     Teuchos::RCP<TreeVector>& soln,
                                     Teuchos::RCP<TreeVector>& soln_dot) {
    (*soln)[0] = S->get_field("pressure", "flow");
    (*soln)[1] = S->get_field("pressure_lambda", "flow");
    (*soln_dot)[0] = S->get_field("pressure_dot", "flow");
    (*soln_dot)[1] = S->get_field("pressure_lambda_dot", "flow");
  };

  /* ******************************************************************
   *  Pointer copy of solution to state
   *********************************************************************/
  void RichardsPK::solution_to_state(Teuchos::RCP<TreeVector>& soln,
                                     Teuchos::RCP<State>& S) {
    S->set_field_pointer("pressure", "flow", (*soln)[0]);
    S->set_field_pointer("pressure_lambda", "flow", (*soln)[1]);
  };

  /* ******************************************************************
   *  Pointer copy of solution to state
   *********************************************************************/
  void RichardsPK::solution_to_state(Teuchos::RCP<TreeVector>& soln,
                                     Teuchos::RCP<TreeVector>& soln_dot,
                                     Teuchos::RCP<State>& S) {
    S->set_field_pointer("pressure", "flow", (*soln)[0]);
    S->set_field_pointer("pressure_lambda", "flow", (*soln)[1]);
    S->set_field_pointer("pressure_dot", "flow", (*soln_dot)[0]);
    S->set_field_pointer("pressure_lambda_dot", "flow", (*soln_dot)[1]);
  };

  /* ******************************************************************
   *  advance transient by a step size dt
   *********************************************************************/
  bool RichardsPK::advance(double dt) {
    h_ = dt;
    state_to_solution(S_next_, solution_);
    time_stepper_->bdf2_step(h_, 0.0, solution_, hnext_);
    time_stepper_->commit_solution(h_, solution_);
    time_stepper_->write_bdf2_stepping_statistics();

    // QUESTION FOR MARKUS -- is it safe to assume that ALL time integration
    // schemes will call fun(solution_) with the final solution?  If not, we
    // must update saturation and flux manually after the solution is
    // acquired.

    return false; // bdf2_step subcycles, so we have always succeeded
  };


  /* ******************************************************************
   *  update diagnostic variables for an I/O call
   *********************************************************************/
  void RichardsPK::calculate_diagnostics(Teuchos::RCP<State>& S) {
    Teuchos::RCP<Epetra_MultiVector> darcy_flux = S->get_field("darcy_flux", "flow");
    Teuchos::RCP<Epetra_MultiVector> darcy_velocity = S->get_field("darcy_velocity", "flow");
    problem_->DeriveDarcyVelocity(*darcy_flux, darcy_velocity);
  };

  ////// STOPPED HERE ////////


  // BDF2 interface
  /* ******************************************************************
   *  computes the non-linear functional f = f(t,u,udot)
   *********************************************************************/
  void RichardsPK::fun(double t, Teuchos::RCP<const TreeVector> soln,
                       Teuchos::RCP<const TreeVector> soln_dot,
                       Teuchos::RCP<TreeVector> f) {
  };

  /* ******************************************************************
   * applies preconditioner to u and returns the result in Pu
   *********************************************************************/
  void RichardsPK::precon(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<TreeVector> Pu);

  /* ******************************************************************
   * computes a norm on u-du and returns the result
   *********************************************************************/
  double RichardsPK::enorm(Teuchos::RCP<const TreeVector> u, Teuchos::RCP<const TreeVector> du);

  /* ******************************************************************
   * updates the preconditioner
   *********************************************************************/
  void RichardsPK::update_precon(double t, Teuchos::RCP<const TreeVector> up, double h, int* errc);


} // close namespace Amanzi
