#include "SteadyState_Richards_PK.hpp"

#include "RichardsProblem.hpp"
#include "RichardsNoxInterface.hpp"

#include "NOX.H"
#include "NOX_Epetra.H"
#include "NOX_GlobalData.H"
#include "NOX_Direction_UserDefinedFactory.H"

#include "NKA.H"
#include "NKADirFactory.H"
#include "NKADirection.H"

SteadyState_Richards_PK::SteadyState_Richards_PK(Teuchos::ParameterList &plist, const Teuchos::RCP<const Flow_State> FS_) : FS(FS_)
{
  // Create the flow boundary conditions object.
  Teuchos::ParameterList bc_plist = plist.sublist("Flow BC");
  bc = Teuchos::rcp<FlowBC>(new FlowBC(bc_plist, FS->mesh()));

  // Create the Richards flow problem.
  problem = new RichardsProblem(FS->mesh(), plist.sublist("Richards Problem"), bc);

  // Create the solution vectors.
  solution = new Epetra_Vector(problem->Map());
  pressure = problem->CreateCellView(*solution);
  darcy_flux = new Epetra_Vector(problem->FaceMap());

  // Get some solver parameters from the flow parameter list.
  max_itr = plist.get<int>("Max Iterations");
  err_tol = plist.get<double>("Error Tolerance");
  precon_freq = plist.get<int>("Preconditioner Update Frequency");

  // Setup the nonlinear solver.
  std::string nl_solver = plist.get<std::string>("Nonlinear Solver");
  if (nl_solver == "JFNK")
    nox_jfnk_setup(nox_param_p, linsol_param_p);
  else if (nl_solver == "NLK")
    nox_nlk_setup(plist, nox_param_p, linsol_param_p);
  else
    throw std::exception();
};

void SteadyState_Richards_PK::nox_jfnk_setup(Teuchos::RCP<Teuchos::ParameterList> &nox_param_p, Teuchos::RCP<Teuchos::ParameterList> &linsol_param_p) const
{
  nox_param_p = Teuchos::rcp<Teuchos::ParameterList>(new Teuchos::ParameterList);
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
  linsol_param_p = Teuchos::rcpFromRef(linsol_param); // we need this later
  linsol_param.set("Aztec Solver", "GMRES");
  //linsol_param.set("Max Iterations", 100);  // INPUT PARAMETER?
  linsol_param.set("Preconditioner", "User Defined");
  linsol_param.set("Preconditioner Reuse Policy", "Recompute");  // INPUT PARAMETER?

  // Set how accurately we solve the linear Newton update problem.
  newton_param.set("Forcing Term Method", "Type 2");
  linsol_param.set("Tolerance", 1.0e-8);  // INPUT PARAMETER?

  // Set the printing parameters in the "Printing" sublist
  Teuchos::ParameterList &print_param = nox_param.sublist("Printing");
  print_param.set("MyPID", problem->Comm().MyPID());
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
}

void SteadyState_Richards_PK::nox_nlk_setup(Teuchos::ParameterList &plist, Teuchos::RCP<Teuchos::ParameterList> &nox_param_p, Teuchos::RCP<Teuchos::ParameterList> &linsol_param_p) const
{
  nox_param_p = Teuchos::rcp<Teuchos::ParameterList>(new Teuchos::ParameterList);
  Teuchos::ParameterList &nox_param = *nox_param_p.get();

  // Use a line search method...
  nox_param.set("Nonlinear Solver", "Line Search Based");
  // with our user-defined NLK method to get the search direction...
  Teuchos::ParameterList &direct_param = nox_param.sublist("Direction");
  direct_param.set("Method", "User Defined");
  // and taking the full step given by NLK.
  Teuchos::ParameterList &search_param = nox_param.sublist("Line Search");
  search_param.set("Method", "Full Step");

  // NKA parameters.
  Teuchos::ParameterList nka_sublist = plist.sublist("NKADirection");

  Teuchos::ParameterList &nka_param = direct_param.sublist("NKADirection");
  nka_param.set("maxv", nka_sublist.get<int>("maxv",5) );
  nka_param.set("vtol", nka_sublist.get<double>("vtol",1.0e-3));

  // Wire-in the NKA direction factory.
  Teuchos::RCP<NOX::GlobalData> gd(new NOX::GlobalData(nox_param_p));
  NOX::Epetra::Vector nox_solution(Teuchos::rcpFromRef(*solution), NOX::Epetra::Vector::CreateView);
  Teuchos::RCP<NOX::Direction::UserDefinedFactory> dir_factory(new NKADirFactory(gd, nka_param, nox_solution));
  direct_param.set("User Defined Direction Factory", dir_factory);

  // Set the linear solver parameters; all that matters is the preconditioner.
  Teuchos::ParameterList &linsol_param = nka_param.sublist("Linear Solver");
  linsol_param_p = Teuchos::rcpFromRef(linsol_param); // we need this later
  linsol_param.set("Aztec Solver", "GMRES");
  linsol_param.set("Preconditioner", "User Defined");

  // Set the printing parameters in the "Printing" sublist
  Teuchos::ParameterList &print_param = nox_param.sublist("Printing");
  print_param.set("MyPID", problem->Comm().MyPID());
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
}


SteadyState_Richards_PK::~SteadyState_Richards_PK()
{
  delete darcy_flux;
  delete pressure;
  delete solution;
  delete problem;
};


int SteadyState_Richards_PK::advance_to_steady_state()
{
  // Set problem parameters.
  problem->SetFluidDensity(FS->fluid_density());
  problem->SetFluidViscosity(FS->fluid_viscosity());
  problem->SetPermeability(FS->permeability());
  problem->SetGravity(FS->gravity());

  // Get references to the bits needed to build the NOX linear system.
  Teuchos::RCP<RichardsNoxInterface> nox_interface(new RichardsNoxInterface(problem));
  Teuchos::RCP<NOX::Epetra::Interface::Required> Ireq = nox_interface;
  Teuchos::RCP<NOX::Epetra::Interface::Preconditioner> Iprec = nox_interface;
  Teuchos::RCP<Epetra_Operator> precon = Teuchos::rcpFromRef(problem->Precon());
  NOX::Epetra::Vector nox_solution(Teuchos::rcpFromRef(*solution), NOX::Epetra::Vector::CreateView);
  Teuchos::ParameterList& print_param = nox_param_p->sublist("Printing");

  // Create the NOX "linear system".
  Teuchos::RCP<NOX::Epetra::LinearSystemAztecOO>
      nox_ls(new NOX::Epetra::LinearSystemAztecOO(print_param, *linsol_param_p, Ireq, Iprec, precon, nox_solution));

  // Create the NOX "group".
  Teuchos::RCP<NOX::Epetra::Group>
      group(new NOX::Epetra::Group(print_param, Ireq, nox_solution, nox_ls));

  // Create the NOX convergence criteron.
  Teuchos::RCP<NOX::StatusTest::NormF> absres(new NOX::StatusTest::NormF(err_tol));
  Teuchos::RCP<NOX::StatusTest::MaxIters> maxitr(new NOX::StatusTest::MaxIters(max_itr));
  Teuchos::RCP<NOX::StatusTest::Combo> conv_test(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
  conv_test->addStatusTest(absres);
  conv_test->addStatusTest(maxitr);

  // Finally create the NOX solver.
  Teuchos::RCP<NOX::Solver::Generic>
      solver = NOX::Solver::buildSolver(group, conv_test, nox_param_p);

  // Set the preconditioner update policy.
  nox_interface->setPrecLag(precon_freq);
  nox_interface->resetPrecLagCounter();

  // Solve the system.
  NOX::StatusTest::StatusType status = solver->solve();

  // Extract the solution from the solver.
  const NOX::Epetra::Group& final_group = dynamic_cast<const NOX::Epetra::Group&>(solver->getSolutionGroup());
  const Epetra_Vector& final_soln = (dynamic_cast<const NOX::Epetra::Vector&>(final_group.getX())).getEpetraVector();
  *solution = final_soln;

  // Derive the Darcy fluxes on faces
  double l1_error;
  problem->DeriveDarcyFlux(*solution, *darcy_flux, l1_error);
  std::cout << "L1 norm of the Darcy flux discrepancy = " << l1_error << std::endl;

  if (status == NOX::StatusTest::Converged)
    return 0;
  else
    return 1;
}

void SteadyState_Richards_PK::GetSaturation(Epetra_Vector &s) const
{
  //for (int i = 0; i < s.MyLength(); ++i) s[i] = 1.0;

  problem->DeriveVanGenuchtenSaturation(*pressure, s);

}
