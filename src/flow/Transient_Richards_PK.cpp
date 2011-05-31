#include "Transient_Richards_PK.hpp"

#include "RichardsProblem.hpp"
#include "RichardsNoxInterface.hpp"

#include "NOX.H"
#include "NOX_Epetra.H"
#include "NOX_GlobalData.H"
#include "NOX_Direction_UserDefinedFactory.H"

#include "NKA.H"
#include "NKADirFactory.H"
#include "NKADirection.H"

Transient_Richards_PK::Transient_Richards_PK(Teuchos::ParameterList &plist, const Teuchos::RCP<const Flow_State> FS_) : FS(FS_)
{
  // Create the flow boundary conditions object.
  Teuchos::ParameterList bc_plist = plist.sublist("Flow BC");
  bc = Teuchos::rcp<FlowBC>(new FlowBC(bc_plist, FS->mesh()));

  // Create the Richards flow problem.
  problem = new RichardsProblem(FS->mesh(), plist.sublist("Richards Problem"), bc);

  // Create the solution vectors.
  solution = new Epetra_Vector(problem->Map());
  pressure = problem->CreateCellView(*solution);
  richards_flux = new Epetra_Vector(problem->FaceMap());

  // Get some solver parameters from the flow parameter list.
  max_itr = plist.get<int>("Max Iterations");
  err_tol = plist.get<double>("Error Tolerance");
  precon_freq = plist.get<int>("Preconditioner Update Frequency");

};

Transient_Richards_PK::~Transient_Richards_PK()
{
  delete richards_flux;
  delete pressure;
  delete solution;
  delete problem;
};


int Transient_Richards_PK::advance()
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

  // Derive the Richards fluxes on faces
  double l1_error;
  problem->DeriveDarcyFlux(*solution, *richards_flux, l1_error);
  std::cout << "L1 norm of the Richards flux discrepancy = " << l1_error << std::endl;

  if (status == NOX::StatusTest::Converged)
    return 0;
  else
    return 1;
}

void Transient_Richards_PK::GetSaturation(Epetra_Vector &s) const
{
  //for (int i = 0; i < s.MyLength(); ++i) s[i] = 1.0;

  problem->DeriveVanGenuchtenSaturation(*pressure, s);

}
