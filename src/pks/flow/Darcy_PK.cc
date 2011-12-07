/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

#include "Darcy_PK.hh"
#include "Mesh.hh"
#include "DarcyProblem.hpp"

namespace Amanzi {

Darcy_PK::Darcy_PK(Teuchos::ParameterList &list, Teuchos::RCP<State> S) {
  // Create the flow boundary conditions object.
  Teuchos::ParameterList bc_param_list = list.sublist("Flow BC");
  bc = Teuchos::rcp<FlowBC>(new FlowBC(bc_param_list, S->get_mesh_maps()));

  // Create the Darcy flow problem.
  Teuchos::ParameterList darcy_plist = list.sublist("Darcy Problem");
  problem = new DarcyProblem(S->get_mesh_maps(),darcy_plist, bc);

  // Require fields in the state.
  S->require_field("pressure", Amanzi::AmanziMesh::CELL, "flow");
  S->require_field("darcy_flux", Amanzi::AmanziMesh::FACE, "flow");
  S->require_field("permeability", Amanzi::AmanziMesh::CELL);

  // Create the solution vectors.
  solution = new Epetra_Vector(problem->Map());
  pressure = problem->CreateCellView(*solution);
  darcy_flux = new Epetra_Vector(problem->FaceMap());

  // Create the linear solver and do the preliminary wiring.
  solver = new AztecOO;
  solver->SetUserOperator(&(problem->Matvec()));
  solver->SetPrecOperator(&(problem->Precon()));
  solver->SetAztecOption(AZ_solver, AZ_cg);

  // Get some solver parameters from the flow parameter list.
  max_itr = list.get<int>("Max Iterations");
  err_tol = list.get<double>("Error Tolerance");
};

Darcy_PK::~Darcy_PK() {
  delete problem;
  delete solver;
  delete darcy_flux;
  delete pressure;
  delete solution;
};

int Darcy_PK::advance_to_steady_state(Teuchos::RCP<State> S) {
  // Set problem parameters.
  problem->SetFluidDensity(*S->get_density());
  problem->SetFluidViscosity(*S->get_viscosity());
  problem->SetGravity(*S->get_gravity());
  problem->SetPermeability(*(*S->get_field("permeability"))(0));

  // Perform the final "assembly" of the problem.
  problem->Assemble();

  // Define the RHS for the system.
  Epetra_Vector b(problem->RHS()); // make a copy, Aztec00 will muck with it
  solver->SetRHS(&b);

  // Define the initial solution guess; overwritten with solution.
  solver->SetLHS(solution);

  // Solve the system.
  solver->Iterate(max_itr, err_tol);
  std::cout << "Darcy solver performed " << solver->NumIters() << " iterations." << std::endl
            << "Norm of true residual = " << solver->TrueResidual() << std::endl;

  // Derive the Darcy fluxes on faces
  double l1_error;
  problem->DeriveDarcyFlux(*solution, *darcy_flux, l1_error);
  std::cout << "L1 norm of the Darcy flux discrepancy = " << l1_error << std::endl;

  // Update the state with these values
  S->set_field("pressure", "flow", *pressure);
  S->set_field("darcy_flux", "flow", *darcy_flux);
  return 0;
}

void Darcy_PK::initialize_state(Teuchos::RCP<State> S) {
  // Calculate the steady state solution.
  int zero = advance_to_steady_state(S);
}

} // close namespace Amanzi
