/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

#include "Permafrost_PK.hh"
#include "PermafrostModelEvaluator.hh"

namespace Amanzi {

  Permafrost_PK::Permafrost_PK(Teuchos::ParameterList &plist_, Teuchos::RCP<State> &S) :
    plist(plist_) {
  // Require fields in the state
  // primary variables: pressure
  S->require_field("pressure", Amanzi::AmanziMesh::CELL, "flow");
  S->get_field_record("pressure")->set_io_vis(true);

  // technically primary variable, is MFD considered a first order form?
  S->require_field("darcy_flux", Amanzi::AmanziMesh::FACE, "flow");

  // secondary variables, saturations
  S->require_field("water saturation", Amanzi::AmanziMesh::CELL, "flow");
  S->get_field_record("water saturation")->set_io_vis(true);
  S->require_field("gas saturation", Amanzi::AmanziMesh::CELL, "flow");
  S->get_field_record("gas saturation")->set_io_vis(true);
  S->require_field("ice saturation", Amanzi::AmanziMesh::CELL, "flow");
  S->get_field_record("ice saturation")->set_io_vis(true);

  // required independent variables, not owned by flow
  // temperature should be owned by a thermal PK, or initialized by state
  S->require_field("temperature", Amanzi::AmanziMesh::CELL);

  // permeability likely owned by no one
  S->require_field("permeability", Amanzi::AmanziMesh::CELL);

  // porosity owned by no one or by a subsidence PK?  The point is
  // we don't care!
  S->require_field("porosity", Amanzi::AmanziMesh::CELL);

  // diagnostics/io
  S->require_field("darcy_velocity", Amanzi::AmanziMesh::CELL, "flow", 3);
  S->get_field_record("darcy_velocity")->set_io_vis(true);
  std::vector<std::string> subfieldnames;
  subfieldnames.push_back("Liquid X-Velocity");
  subfieldnames.push_back("Liquid Y-Velocity");
  subfieldnames.push_back("Liquid Z-Velocity");
  S->set_subfield_names("darcy_velocity", subfieldnames);

  // create the problem
  Teuchos::ParameterList permafrost_plist = plist.sublist("Permafrost Problem");
  problem = Teuchos::rcp(new PermafrostProblem(S->get_mesh_maps(), permafrost_plist));

  // Create the solution vectors.
  solution = Teuchos::rcp(new Epetra_Vector(problem->Map()));
  pressure_cells = Teuchos::rcp(problem->CreateCellView(*solution));
  pressure_faces = Teuchos::rcp(problem->CreateFaceView(*solution));
  darcy_flux = Teuchos::rcp(new Epetra_Vector(problem->FaceMap()));

  // and the model evaluator
  Teuchos::ParameterList &rme_list = permafrost_plist.sublist("Permafrost model evaluator");
  RME = Teuchos::rcp(new PermafrostModelEvaluator(problem, rme_list));

  // then the BDF2 solver
  Teuchos::RCP<Teuchos::ParameterList> bdf2_list_p(new Teuchos::ParameterList(permafrost_plist.sublist("Time integrator")));
  time_stepper = Teuchos::rcp(new BDF2::Dae(*RME, problem->Map()));
  time_stepper->setParameterList(bdf2_list_p);
};

void Permafrost_PK::initialize(Teuchos::RCP<State> &S) {
  Teuchos::ParameterList permafrost_plist = plist.sublist("Permafrost Problem");

  steady_state = permafrost_plist.get<bool>("Steady state", false);
  if (steady_state) {
    ss_t0 = permafrost_plist.get<double>("Steady state calculation initial time");
    ss_t1 = permafrost_plist.get<double>("Steady state calculation final time");
    h0 = permafrost_plist.get<double>("Steady state calculation initial time step");
    height0 =  permafrost_plist.get<double>("Steady state calculation initial hydrostatic pressure height", 0.0);
  } else {
    h0 = permafrost_plist.get<double>("Initial time step");
    height0 =  permafrost_plist.get<double>("Initial hydrostatic pressure height", 0.0);
  }
  // set the intial timestep size
  h = h0;
  hnext = h0;

  // Set the problem independent variables from state
  problem->SetFluidDensity(*S->get_density());
  problem->SetFluidViscosity(*S->get_viscosity());
  problem->SetGravity(*S->get_gravity());

  problem->SetPermeability(*(*S->get_field("permeability"))(0));
  problem->SetPorosity(*(*S->get_field("porosity"))(0));

  // initialize problem, including data needed to construct
  // BCs, water models, etc.
  problem->InitializeProblem(permafrost_plist);

  // IC for pressure
  problem->SetInitialPressureProfileCells(height0,pressure_cells);
  problem->SetInitialPressureProfileFaces(height0,pressure_faces);
  S->set_field("pressure", "flow", *pressure_cells);

  // IC for saturation -- require it to be set by the state

  // IC for flux
  double l1_error;
  problem->DeriveDarcyFlux(*solution, *darcy_flux, l1_error);
  S->set_field("darcy_flux", "flow", *darcy_flux);

  // initialize timestepping
  if (steady_state) {
    S->set_time(ss_t0);
  }

  double t0 = S->get_time();
  Epetra_Vector udot(problem->Map());
  problem->Compute_udot(t0, *solution, udot);
  time_stepper->set_initial_state(t0, *solution, udot);

  // set up preconditioner
  int errc;
  RME->update_precon(t0, *solution, h, errc);

  // if steady state solve
  if (steady_state) {
    // NOTE: this does NOT get (nor need) state, to gaurantee it cannot alter State
    advance_to_steady_state();

    // update state
    S->set_field("pressure", "flow", *pressure_cells);
    S->set_field("darcy_flux", "flow", *darcy_flux);
    S->set_time(ss_t1);
  }

  // IC for velocity
  Teuchos::RCP<Epetra_MultiVector> velocity = S->get_field("darcy_velocity", "flow");
  problem->DeriveDarcyVelocity(*solution, *velocity);
};

void Permafrost_PK::advance_to_steady_state()
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

  // Derive the Permafrost fluxes on faces
  double l1_error;
  problem->DeriveDarcyFlux(*solution, *darcy_flux, l1_error);
  std::cout << "L1 norm of the Permafrost flux discrepancy = " << l1_error << std::endl;
}

// NOTE: this should get split into three parts --
//    1. set problem data from state S0
//    2. advance transient (which in no way can alter either S0 or S1)
//    3. set state S1 with solution from problem
// much like advance_to_steady_state() works
bool Permafrost_PK::advance_transient(double dt, const Teuchos::RCP<State> &S0,
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
  problem->DeriveDarcyFlux(*solution, *darcy_flux, l1_error);
  std::cout << "L1 norm of the Permafrost flux discrepancy = " << l1_error << std::endl;
  S1->set_field("darcy_flux", "flow", *darcy_flux);

  // update velocity DOES THIS NEED TO BE DONE NOW?  save for commit_state?
  Teuchos::RCP<Epetra_MultiVector> velocity = S1->get_field("darcy_velocity", "flow");
  problem->DeriveDarcyVelocity(*solution, *velocity);

  // update saturation
  Teuchos::RCP<Epetra_MultiVector> saturation = S1->get_field("water saturation", "flow");
  problem->DeriveSaturation(*pressure_cells, *(*saturation)(0));

  return false; // bdf2_step subcycles, so we have always succeeded
};
} // close namespace Amanzi
