#include "Transient_Richards_PK.hpp"

#include "RichardsProblem.hpp"

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

  // create the time stepper...

  // first the Richards model evaluator
  Teuchos::ParameterList &rme_list = plist.sublist("Richards model evaluator");
  RME = new RichardsModelEvaluator(problem, rme_list, problem->Map(), FS);  

  // then the BDF2 solver
  Teuchos::RCP<Teuchos::ParameterList> bdf2_list_p(new Teuchos::ParameterList(plist.sublist("Time integrator")));

  time_stepper = new BDF2::Dae(*RME, problem->Map());
  time_stepper->setParameterList(bdf2_list_p);

};

Transient_Richards_PK::~Transient_Richards_PK()
{
  delete richards_flux;
  delete pressure;
  delete solution;
  delete problem;
};


int Transient_Richards_PK::advance_to_steady_state()
{
  // Set problem parameters.
  problem->SetFluidDensity(FS->fluid_density());
  problem->SetFluidViscosity(FS->fluid_viscosity());
  problem->SetPermeability(FS->permeability());
  problem->SetGravity(FS->gravity());
  problem->SetFlowState(FS);

  double t0 = 0.0;
  double t1 = 100.0;
  double h = 1e-5;
  double hnext;


  // create udot
  Epetra_Vector udot(problem->Map());
  problem->Compute_udot(t0,  *solution, udot);

  time_stepper->set_initial_state(t0, *solution, udot);

  int errc;
  RME->update_precon(t0, *solution, h, errc);

  // iterate
  int i = 0;
  double tlast = t0;

  
  do {
    
    time_stepper->bdf2_step(h,0.0,20,*solution,hnext);
    
    time_stepper->commit_solution(h,*solution);

    // time_stepper->advance_time(h);

    // update the state, but only the cell values of pressure
    // FS->update_pressure( * problem->CreateCellView(*solution) );

    time_stepper->write_bdf2_stepping_statistics();

    h = hnext;
    i++;

    tlast=time_stepper->most_recent_time();


  } while (t1 >= tlast);    
  


  // Derive the Richards fluxes on faces
  double l1_error;
  problem->DeriveDarcyFlux(*solution, *richards_flux, l1_error);
  std::cout << "L1 norm of the Richards flux discrepancy = " << l1_error << std::endl;

}

void Transient_Richards_PK::GetSaturation(Epetra_Vector &s) const
{
  //for (int i = 0; i < s.MyLength(); ++i) s[i] = 1.0;

  problem->DeriveVanGenuchtenSaturation(*pressure, s);

}
