/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

#include "Darcy_PK.hh"
#include "DarcyProblem.hpp"

namespace Amanzi {

Darcy_PK::Darcy_PK(Teuchos::ParameterList &plist_, Teuchos::RCP<State> &S) :
    plist(plist_) {

  // Require fields in the state.
  S->require_field("pressure", Amanzi::AmanziMesh::CELL, "flow");
  S->get_field_record("pressure")->set_io_vis(true);

  S->require_field("darcy_flux", Amanzi::AmanziMesh::FACE, "flow");

  S->require_field("darcy_velocity", Amanzi::AmanziMesh::CELL, "flow", 3);
  S->get_field_record("darcy_velocity")->set_io_vis(true);
  std::vector<std::string> subfieldnames;
  subfieldnames.push_back("Liquid X-Velocity");
  subfieldnames.push_back("Liquid Y-Velocity");
  subfieldnames.push_back("Liquid Z-Velocity");
  S->set_subfield_names("darcy_velocity", subfieldnames);

  S->require_field("permeability", Amanzi::AmanziMesh::CELL);

  // Create the Darcy flow problem.
  Teuchos::ParameterList darcy_plist = plist.sublist("Darcy Problem");
  problem = Teuchos::rcp(new DarcyProblem(S->get_mesh_maps(), darcy_plist));

  // Create the solution vectors.
  solution = Teuchos::rcp(new Epetra_Vector(problem->Map()));
  pressure = Teuchos::rcp(problem->CreateCellView(*solution));
  darcy_flux = Teuchos::rcp(new Epetra_Vector(problem->FaceMap()));

  // Create the linear solver
  solver = Teuchos::rcp(new AztecOO);
};

void Darcy_PK::advance_to_steady_state() {

  // Perform the final "assembly" of the problem.
  problem->Assemble();

  // Define the RHS for the system.
  Epetra_Vector b(problem->RHS()); // make a copy, Aztec00 will muck with it
  solver->SetRHS(&b);

  // Define the initial solution guess; overwritten with solution.
  solver->SetLHS(&(*solution));

  // Solve the system.
  solver->Iterate(max_itr, err_tol);
  std::cout << "Darcy solver performed " << solver->NumIters() << " iterations." << std::endl
            << "Norm of true residual = " << solver->TrueResidual() << std::endl;

  // Derive the Darcy fluxes on faces
  double l1_error;
  problem->DeriveDarcyFlux(*solution, *darcy_flux, l1_error);
  std::cout << "L1 norm of the Darcy flux discrepancy = " << l1_error << std::endl;
};

void Darcy_PK::initialize(Teuchos::RCP<State> &S) {
  // Get some solver parameters from the flow parameter list.
  max_itr = plist.get<int>("Max Iterations");
  err_tol = plist.get<double>("Error Tolerance");

  // Set problem parameters.
  problem->SetFluidDensity(*S->get_density());
  problem->SetFluidViscosity(*S->get_viscosity());
  problem->SetGravity(*S->get_gravity());
  problem->SetPermeability(*(*S->get_field("permeability"))(0));

  // initialize problem, including data needed to construct
  // BCs, water models, etc.
  Teuchos::ParameterList darcy_plist = plist.sublist("Darcy Problem");
  problem->InitializeProblem(darcy_plist);

  // Initialize the solver
  solver->SetUserOperator(&(problem->Matvec()));
  solver->SetPrecOperator(&(problem->Precon()));
  solver->SetAztecOption(AZ_solver, AZ_cg);

  // Calculate the steady state solution.
  // NOTE: this does NOT get (nor need) state, to gaurantee it cannot alter State
  advance_to_steady_state();

  // Update the state with these values
  S->set_field("pressure", "flow", *pressure);
  S->set_field("darcy_flux", "flow", *darcy_flux);
  Teuchos::RCP<Epetra_MultiVector> velocity = S->get_field("darcy_velocity", "flow");
  problem->DeriveDarcyVelocity(*solution, *velocity);
};

} // close namespace Amanzi
