/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

#include "DarcyPK.hh"
#include "DarcyProblem.hpp"

namespace Amanzi {

DarcyPK::DarcyPK(Teuchos::ParameterList &flow_plist, Teuchos::RCP<State> &S,
                 Teuchos::RCP<TreeVector>& soln) :
    flow_plist(plist_) {

  // Require fields in the state.
  S->require_field("pressure", FIELD_LOCATION_CELL, "flow");
  S->get_field_record("pressure")->set_io_vis(true);
  Teuchos::RCP<Epetra_MultiVector> soln_prs =
    Teuchos::rcp(new Epetra_MultiVector(*S->get_field("pressure")));
  soln->PushBack(soln_prs);

  S->require_field("pressure_lambda", FIELD_LOCATION_FACE, "flow");
  Teuchos::RCP<Epetra_MultiVector> lambda_soln_prs =
    Teuchos::rcp(new Epetra_MultiVector(*S->get_field("pressure_lambda")));
  soln->PushBack(lambda_soln_prs);

  // Require flux for transport/etc.
  S->require_field("darcy_flux", FIELD_LOCATION_FACE, "flow");

  // Diagnostic fields.
  S->require_field("darcy_velocity", FIELD_LOCATION_CELL, "flow", 3);
  S->get_field_record("darcy_velocity")->set_io_vis(true);
  std::vector<std::string> subfieldnames;
  subfieldnames.push_back("Liquid X-Velocity");
  subfieldnames.push_back("Liquid Y-Velocity");
  subfieldnames.push_back("Liquid Z-Velocity");
  S->set_subfield_names("darcy_velocity", subfieldnames);

  // required constants/independent variables
  S->require_field("permeability", FIELD_LOCATION_CELL);

  // Create the Darcy flow problem.
  Teuchos::ParameterList darcy_plist = plist.sublist("Darcy Problem");
  problem_ = Teuchos::rcp(new DarcyProblem(S->get_mesh_maps(), darcy_plist));

  // Create the linear solver
  solver_ = Teuchos::rcp(new AztecOO);
};

// Private method to calculate the Darcy solution.
void DarcyPK::advance_to_steady_state(Teuchos::RCP<TreeVector>& soln) {
  // Perform the final "assembly" of the problem.
  problem_->Assemble();

  // Define the RHS for the system.
  Epetra_Vector b(problem_->RHS()); // make a copy, Aztec00 will muck with it
  solver_->SetRHS(&b);

  // Define the initial solution guess; overwritten with solution.
  Teuchos::RCP<TreeVector> pressure;
  int subvec_not_found = soln->SubVector("pressure", pressure);
  if (subvec_not_found) {
    Errors::Message message("MPC: vector structure does not match PK structure");
    Exceptions::amanzi_throw(message);
  }
  solver_->SetLHS(&(*pressure)); // this is bad... passing raw pointer...

  // Solve the system.
  solver_->Iterate(max_itr_, err_tol_);
  std::cout << "Darcy solver performed " << solver_->NumIters() << " iterations."
            << std::endl << "Norm of true residual = " << solver_->TrueResidual()
            << std::endl;
};

void DarcyPK::initialize(Teuchos::RCP<State> &S, Teuchos::RCP<TreeVector>& soln) {
  // Get some solver parameters from the flow parameter list.
  max_itr_ = flow_plist_.get<int>("Max Iterations");
  err_tol_ = flow_plist_.get<double>("Error Tolerance");

  // Set problem parameters.
  problem_->SetFluidDensity(*S->get_density());
  problem_->SetFluidViscosity(*S->get_viscosity());
  problem_->SetGravity(*S->get_gravity());
  problem_->SetPermeability(*(*S->get_field("permeability"))(0));

  // initialize problem, including data needed to construct
  // BCs, water models, etc.
  Teuchos::ParameterList darcy_plist = flow_plist_.sublist("Darcy Problem");
  problem_->InitializeProblem(darcy_plist);

  // Initialize the solver
  solver_->SetUserOperator(&(problem_->Matvec()));
  solver_->SetPrecOperator(&(problem_->Precon()));
  solver_->SetAztecOption(AZ_solver, AZ_cg);

  // Calculate the steady state solution.
  // NOTE: this does NOT get (nor need) state, to gaurantee it cannot alter State
  advance_to_steady_state();


  // Derive the Darcy fluxes on faces
  double l1_error;
  problem_->DeriveDarcyFlux(*pressure, *darcy_flux, l1_error);
  std::cout << "L1 norm of the Darcy flux discrepancy = " << l1_error << std::endl;

  // Update the state with these values
  S->set_field("pressure", "flow", *pressure);
  S->set_field("darcy_flux", "flow", *darcy_flux);
  Teuchos::RCP<Epetra_MultiVector> velocity = S->get_field("darcy_velocity", "flow");
  problem_->DeriveDarcyVelocity(*solution, *velocity);
};

} // close namespace Amanzi
