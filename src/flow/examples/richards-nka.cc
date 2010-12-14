#include "Mesh_maps_moab.hh"

#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Epetra_MpiComm.h"
#include "mpi.h"

#include "NOX.H"
#include "NOX_Epetra.H"
#include "NOX_GlobalData.H"
#include "NOX_Direction_UserDefinedFactory.H"

#include "NKA.H"
#include "NKADirFactory.H"
#include "NKADirection.H"

#include "RichardsProblem.hpp"
#include "gmv_mesh.hh"

enum Side { LEFT, RIGHT, FRONT, BACK, BOTTOM, TOP };

int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  Epetra_MpiComm comm(MPI_COMM_WORLD);
  
  // MESH
  Teuchos::RCP<Mesh_maps_moab> mesh(new Mesh_maps_moab(argv[1], MPI_COMM_WORLD));
  
  // BOUNDARY CONDITIONS
  Teuchos::ParameterList bc_params;
  bc_params.set("number of BCs", 6);
  Teuchos::ParameterList &bc_left = bc_params.sublist("BC00");
    bc_left.set("Side set ID", 4);
    bc_left.set("Type", "Pressure Constant");
    bc_left.set("BC value", 1.0);
  Teuchos::ParameterList &bc_right = bc_params.sublist("BC01");
    bc_right.set("Side set ID", 2);
    bc_right.set("Type", "Pressure Constant");
    bc_right.set("BC value", 0.0);
  Teuchos::ParameterList &bc_front = bc_params.sublist("BC02");
    bc_front.set("Side set ID", 1);
    bc_front.set("Type", "No Flow");
  Teuchos::ParameterList &bc_back = bc_params.sublist("BC03");
    bc_back.set("Side set ID", 3);
    bc_back.set("Type", "No Flow");
  Teuchos::ParameterList &bc_bottom = bc_params.sublist("BC04");
    bc_bottom.set("Side set ID", 5);
    bc_bottom.set("Type", "No Flow");
  Teuchos::ParameterList &bc_top = bc_params.sublist("BC05");
    bc_top.set("Side set ID", 6);
    bc_top.set("Type", "No Flow");
  Teuchos::RCP<FlowBC> bc(new FlowBC(bc_params, mesh));
  
  // PROBLEM
  Teuchos::ParameterList pl;
  RichardsProblem problem(mesh, pl, bc);

  // MODEL PARAMETERS
  problem.SetFluidDensity(1.0);
  problem.SetFluidViscosity(1.0);
  problem.SetPermeability(1.0);
  double g[3] = { 0.0 };
  problem.SetGravity(g);
  
  // Create the NOX parameter list.
  Teuchos::RCP<Teuchos::ParameterList> nox_param_p(new Teuchos::ParameterList);
  Teuchos::ParameterList &nox_param = *nox_param_p.get();
 
  // we need a GlobalData variable fo initialize the 
  // NKA Direction factory later on
  Teuchos::RCP<NOX::GlobalData> gd = Teuchos::rcp(new NOX::GlobalData(nox_param_p));


  // Use a line search method...
  nox_param.set("Nonlinear Solver", "Line Search Based");
  
  // with get the search direction...
  Teuchos::ParameterList &direct_param = nox_param.sublist("Direction");
  direct_param.set("Method", "User Defined");
  
  // these are parameters that the NKA direction class
  // will look at
  Teuchos::ParameterList &mydirParams = nox_param.sublist("NKADirection");
  mydirParams.set("maxv", 4);
  mydirParams.set("vtol", 5e-2);
  

  // and taking the full step given by the nonlinear solver
  Teuchos::ParameterList &search_param = nox_param.sublist("Line Search");
  search_param.set("Method", "Full Step");
  
  Teuchos::ParameterList &linsol_param = direct_param.sublist("Linear Solver");
  linsol_param.set("Preconditioner", "User Defined");

  // Set the printing parameters in the "Printing" sublist
  Teuchos::ParameterList &print_param = nox_param.sublist("Printing");
  print_param.set("MyPID", comm.MyPID()); 
  print_param.set("Output Precision", 5);
  print_param.set("Output Processor", 0);
  print_param.set("Output Information", 
			   NOX::Utils::OuterIteration + 
			   NOX::Utils::OuterIterationStatusTest + 
			   NOX::Utils::InnerIteration +
			   NOX::Utils::LinearSolverDetails +
			   NOX::Utils::Parameters + 
			   NOX::Utils::Details + 
			   NOX::Utils::Warning +
			   NOX::Utils::Error +
                           NOX::Utils::Debug);


  Teuchos::RCP<NOX::Epetra::Interface::Required> Ireq = Teuchos::rcpFromRef(problem.NoxReq());                           
  Teuchos::RCP<NOX::Epetra::Interface::Preconditioner> Iprec = Teuchos::rcpFromRef(problem.NoxPrecon());
  Teuchos::RCP<Epetra_Operator> precon = Teuchos::rcpFromRef(problem.Precon());
  
  //Teuchos::RCP<Epetra_Vector> solution(new Epetra_Vector(problem.Map()));
  Epetra_Vector solution(problem.Map());
  NOX::Epetra::Vector nox_solution(Teuchos::rcpFromRef(solution), NOX::Epetra::Vector::CreateView);
  
  // I don't think this should be necessary; NOX should do this at the start of the iteration.
  //problem.ComputePrecon(solution);
  
  // Create the NOX "linear system".
  Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO>
      nox_ls(new NOX::Epetra::LinearSystemAztecOO(print_param, linsol_param, Ireq, Iprec, precon, nox_solution));


  // create our NKA direction factory
  Teuchos::RCP<NOX::Direction::UserDefinedFactory> dir_factory 
    = Teuchos::rcp(new NKADirFactory(gd, mydirParams, nox_solution));
  direct_param.set("User Defined Direction Factory", dir_factory);

  // Create the NOX "group".
  Teuchos::RCP<NOX::Epetra::Group> group(new NOX::Epetra::Group(print_param, Ireq, nox_solution, nox_ls));

  // Create the NOX convergence criteron.
  Teuchos::RCP<NOX::StatusTest::NormF> abs_res(new NOX::StatusTest::NormF(1.0e-12));
  Teuchos::RCP<NOX::StatusTest::MaxIters> max_itr(new NOX::StatusTest::MaxIters(100));
  Teuchos::RCP<NOX::StatusTest::Combo> conv_test(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
  conv_test->addStatusTest(abs_res);
  conv_test->addStatusTest(max_itr);
  
  // Finally we create the NOX solver for the problem
  Teuchos::RCP<NOX::Solver::Generic> solver = NOX::Solver::buildSolver(group, conv_test, nox_param_p);
  
  std::cout << "ABOUT TO CALL SOLVE ..." << std::endl;
  
  // Solve the nonlinear system
  NOX::StatusTest::StatusType status = solver->solve();
  
  // Extract the solution from the solver.
  const NOX::Epetra::Group& final_group = dynamic_cast<const NOX::Epetra::Group&>(solver->getSolutionGroup());
  const Epetra_Vector& final_soln = (dynamic_cast<const NOX::Epetra::Vector&>(final_group.getX())).getEpetraVector();
  
  // Write a GMV file with the cell pressures.
  std::string gmv_file("pressure.gmv");
  GMV::open_data_file(*mesh, gmv_file);
  GMV::start_data();
  Epetra_Vector &Pcell = *(problem.CreateCellView(final_soln));
  GMV::write_cell_data(Pcell, std::string("pressure"));
  GMV::close_data_file();
  delete &Pcell;
  
  return 0;
}
