#include "Mesh_MOAB.hh"

#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Epetra_MpiComm.h"
#include "mpi.h"

#include "NOX.H"
#include "NOX_Epetra.H"

#include "RichardsProblem.hpp"
#include "RichardsNoxInterface.hpp"
#include "FlowBC.hpp"
#include "gmv_mesh.hh"

int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  Epetra_MpiComm comm(MPI_COMM_WORLD);

  // MESH
  Teuchos::RCP<Amanzi::AmanziMesh::Mesh_MOAB> mesh(new Amanzi::AmanziMesh::Mesh_MOAB(argv[1], MPI_COMM_WORLD));

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
    Teuchos::RCP<Amanzi::FlowBC> bc(new Amanzi::FlowBC(bc_params, mesh));

  // PROBLEM
  Teuchos::ParameterList pl;
  Amanzi::RichardsProblem problem(mesh, pl, bc);

  // MODEL PARAMETERS
  problem.SetFluidDensity(1.0);
  problem.SetFluidViscosity(1.0);
  problem.SetPermeability(1.0);
  double g[3] = { 0.0 };
  problem.SetGravity(g);

  // Initial solution vector and its NOX view.
  Epetra_Vector solution(problem.Map());
  NOX::Epetra::Vector nox_solution(Teuchos::rcpFromRef(solution), NOX::Epetra::Vector::CreateView);

  // Create the NOX parameter list.
  Teuchos::RCP<Teuchos::ParameterList> nox_param_p(new Teuchos::ParameterList);
  Teuchos::ParameterList &nox_param = *nox_param_p.get();

  // Use a line search method...
  nox_param.set("Nonlinear Solver", "Line Search Based");

  // with Newton to get the search direction...
  Teuchos::ParameterList &direct_param = nox_param.sublist("Direction");
  direct_param.set("Method", "Newton");

  // and taking the full step given by Newton.
  Teuchos::ParameterList &search_param = nox_param.sublist("Line Search");
  search_param.set("Method", "Full Step");

  // Set the Newton and linear solver parameters.
  Teuchos::ParameterList &newton_param = direct_param.sublist("Newton");
  Teuchos::ParameterList &linsol_param = newton_param.sublist("Linear Solver");
  linsol_param.set("Aztec Solver", "GMRES");
  linsol_param.set("Max Iterations", 100);
  linsol_param.set("Preconditioner", "User Defined");

  // Set how accurately we solve the linear Newton update problem.
  newton_param.set("Forcing Term Method", "Constant");
  linsol_param.set("Tolerance", 1.0e-8);

  // Set the printing parameters in the "Printing" sublist.
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
			   NOX::Utils::Error);

  // Collect the bits needed to build the NOX linear system.
  Teuchos::RCP<Amanzi::RichardsNoxInterface> nox_interface(new Amanzi::RichardsNoxInterface(&problem));
  Teuchos::RCP<NOX::Epetra::Interface::Required> Ireq = nox_interface;
  Teuchos::RCP<NOX::Epetra::Interface::Preconditioner> Iprec = nox_interface;
  Teuchos::RCP<Epetra_Operator> precon = Teuchos::rcpFromRef(problem.Precon());

  // Create the NOX "linear system".
  Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO>
      nox_ls(new NOX::Epetra::LinearSystemAztecOO(print_param, linsol_param, Ireq, Iprec, precon, nox_solution));

  // Create the NOX "group".
  Teuchos::RCP<NOX::Epetra::Group>
      group(new NOX::Epetra::Group(print_param, Ireq, nox_solution, nox_ls));

  // Create the NOX convergence criteron.
  Teuchos::RCP<NOX::StatusTest::NormF> abs_res(new NOX::StatusTest::NormF(1.0e-12));
  Teuchos::RCP<NOX::StatusTest::MaxIters> max_itr(new NOX::StatusTest::MaxIters(100));
  Teuchos::RCP<NOX::StatusTest::Combo> conv_test(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
  conv_test->addStatusTest(abs_res);
  conv_test->addStatusTest(max_itr);

  // Finally we create the NOX solver for the problem.
  Teuchos::RCP<NOX::Solver::Generic> solver = NOX::Solver::buildSolver(group, conv_test, nox_param_p);

  // Solve the nonlinear system.
  NOX::StatusTest::StatusType status = solver->solve();

  // Extract the solution from the solver.
  const NOX::Epetra::Group& final_group = dynamic_cast<const NOX::Epetra::Group&>(solver->getSolutionGroup());
  const Epetra_Vector& final_soln = (dynamic_cast<const NOX::Epetra::Vector&>(final_group.getX())).getEpetraVector();

  // Write a GMV file with the cell pressures.
  std::string gmv_file("pressure.gmv");
  Amanzi::GMV::open_data_file(*mesh, gmv_file);
  Amanzi::GMV::start_data();
  Epetra_Vector &Pcell = *(problem.CreateCellView(final_soln));
  Amanzi::GMV::write_cell_data(Pcell, std::string("pressure"));
  Amanzi::GMV::close_data_file();
  delete &Pcell;

  return 0;
}
