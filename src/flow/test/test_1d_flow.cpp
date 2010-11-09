//#include "mpi.h"
#include "UnitTest++.h"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "Epetra_SerialComm.h"
#include "Epetra_MpiComm.h"

#include "AztecOO.h"

#include "Mesh_maps_simple.hh"
#include "DiffusionMatrix.hpp"
#include "DiffusionPrecon.hpp"
#include "DarcyProblem.hpp"
#include "DarcyMatvec.hpp"

#include "gmv_mesh.hh"

//int main (int argc, char *argv[])
//{
//  MPI_Init(&argc, &argv);
  
TEST(1D_FLOW) {
  
#ifdef HAVE_MPI
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif
  
  
  //
  // CREATE A SIMPLE BRICK MESH
  
  //Teuchos::RCP<Mesh_maps_base> mesh(new Mesh_maps_simple(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 5, 5, 5, &comm));
  //Teuchos::RCP<Mesh_maps_base> mesh(new Mesh_maps_simple(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 100, 1, 1, &comm));
  Teuchos::RCP<Mesh_maps_base> mesh(new Mesh_maps_simple(0.0, 0.0, 0.0, 1.0, 0.1, 0.1, 50, 1, 1, &comm));
  
  
  //
  // CREATE THE BOUNDARY CONDITIONS
  
  // Define the BC parameter list.
  Teuchos::ParameterList bc_list;
  bc_list.set("number of BCs", 6);
  // Left side (ID 3)
  Teuchos::ParameterList &bc_left = bc_list.sublist("BC00");
  bc_left.set("Side set ID", 3);
  bc_left.set("Type", "Pressure Constant");
  bc_left.set("BC value", 0.0);
  // Right side (ID 1)
  Teuchos::ParameterList &bc_right = bc_list.sublist("BC01");
  bc_right.set("Side set ID", 1);
  bc_right.set("Type", "Pressure Constant");
  bc_right.set("BC value", 1.0);
  // Front side (ID 0)
  Teuchos::ParameterList &bc_front = bc_list.sublist("BC02");
  bc_front.set("Side set ID", 0);
  bc_front.set("Type", "No Flow");
  // Back side (ID 2)
  Teuchos::ParameterList &bc_back = bc_list.sublist("BC03");
  bc_back.set("Side set ID", 2);
  bc_back.set("Type", "No Flow");
  // Bottom side (ID 4)
  Teuchos::ParameterList &bc_bottom = bc_list.sublist("BC04");
  bc_bottom.set("Side set ID", 4);
  bc_bottom.set("Type", "No Flow");
  // Bottom side (ID 5)
  Teuchos::ParameterList &bc_top = bc_list.sublist("BC05");
  bc_top.set("Side set ID", 5);
  bc_top.set("Type", "No Flow");
  
  // Create the flow BCs from the parameter list.
  Teuchos::RCP<FlowBC> bc(new FlowBC(bc_list, mesh));

  
  //
  // CREATE THE DARCY FLOW PROBLEM
  
  DarcyProblem prob(mesh, bc);
  
  prob.SetFluidDensity(1.0);
  prob.SetFluidViscosity(1.0);
  prob.SetPermeability(1.0);
  double g[3] = { 0.0 };
  prob.SetGravity(g);
  
  prob.Assemble();
  
  // Acquire the matvec operator for the system.
  Epetra_Operator &Ax = prob.Matvec();
  
  // Acquire the preconditioning operator for the system.
  Epetra_Operator &precon = prob.Precon();
  
  
  //
  // CREATE THE AZTECOO SOLVER FOR THE LINEAR SYSTEM
  
  AztecOO solver; // an empty solver
  
  // Register the matvec operator (Epetra_Operator)
  solver.SetUserOperator(&Ax);
  
  // Register the preconditioning operator (Epetra_Operator)
  solver.SetPrecOperator(&precon);
  
  // Register the RHS.
  Epetra_Vector B(prob.RHS()); // make a copy, Aztec00 wants to muck with it
  solver.SetRHS(&B);
  
  // Register the initial solution guess; will be overwritten with the solution.
  Epetra_Vector X(prob.Map()); // constructor sets values to zero
  solver.SetLHS(&X);
  
  // Set the solver options.
  solver.SetAztecOption(AZ_solver, AZ_cg); // use CG; the system is symmetric
  //solver.SetAztecOption(AZ_precond, AZ_none);
  
  
  //
  // SOLVE THE SYSTEM
  
  solver.Iterate(100, 1.0e-8);
  std::cout << "Solver performed " << solver.NumIters() << " iterations." << std::endl
            << "Norm of true residual = " << solver.TrueResidual() << std::endl;
  
  // Print the solution; contains both cell and face pressures.
  //std::cout << X << std::endl;
  
  Epetra_Vector &Pcell = *(prob.CreateCellView(X));
  //std::cout << Pcell << std::endl;
  
  // Write a GMV file with the cell pressures.
  std::string gmv_file("flow.gmv");
  GMV::open_data_file(*mesh, gmv_file);
  GMV::start_data();
  GMV::write_cell_data(Pcell, std::string("pressure"));
  GMV::close_data_file();
  delete &Pcell;
  
//  MPI_Finalize();
//  return 0;
}
