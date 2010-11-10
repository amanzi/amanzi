#include "Flow_PK.hpp"

#include "DarcyProblem.hpp"

Flow_PK::Flow_PK(Teuchos::ParameterList &list, const Teuchos::RCP<const Flow_State> FS_) : FS(FS_)
{
  // Create the flow boundary conditions object.
  Teuchos::ParameterList bc_param_list = list.sublist("Flow BC");
  bc = Teuchos::rcp<FlowBC>(new FlowBC(bc_param_list, FS->mesh()));
  
  // Create the Darcy flow problem.
  problem = new DarcyProblem(FS->mesh(), bc);
  
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

Flow_PK::~Flow_PK()
{ 
  delete solver;
  delete darcy_flux;
  delete pressure;
  delete problem;
};


int Flow_PK::advance()
{
  // Set problem parameters.
  problem->SetFluidDensity(FS->fluid_density());
  problem->SetFluidViscosity(FS->fluid_viscosity());
  problem->SetPermeability(FS->permeability());
  problem->SetGravity(FS->gravity());
  
  // Temporary ...
  problem->SetFluidDensity(1.0);
  problem->SetFluidViscosity(1.0);
  problem->SetPermeability(1.0);
  double g[3] = {0.0, 0.0, -1.0};
  problem->SetGravity(g);
  
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
