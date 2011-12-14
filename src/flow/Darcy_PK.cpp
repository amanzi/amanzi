#include "Darcy_PK.hpp"

#include "DarcyProblem.hpp"

namespace Amanzi
{

Darcy_PK::Darcy_PK(Teuchos::ParameterList &list, const Teuchos::RCP<const Flow_State> FS_) : FS(FS_)
{
  // Add some parameters to the Darcy problem constructor parameter list.
  Teuchos::ParameterList &dp_list = list.sublist("Darcy Problem");
  dp_list.set("fluid density", FS->fluid_density());
  dp_list.set("fluid viscosity", FS->fluid_viscosity());
  const double *gravity = FS->gravity();
  //TODO: assuming gravity[0] = gravity[1] = 0 -- needs to be reconciled somehow
  dp_list.set("gravity", -gravity[2]);
  
  // Create the Darcy flow problem.
  Teuchos::ParameterList darcy_plist = list.sublist("Darcy Problem");
  problem = new DarcyProblem(FS->mesh(),darcy_plist);

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


Darcy_PK::~Darcy_PK()
{
  delete problem;
  delete solver;
  delete darcy_flux;
  delete pressure; 
  delete solution; 
};


int Darcy_PK::advance_to_steady_state()
{
  // Set problem parameters.
  problem->SetPermeability(FS->permeability());

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

  return 0;
}

} // close namespace Amanzi
