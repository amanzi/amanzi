#include "Transient_Richards_PK.hpp"

#include "RichardsProblem.hpp"

namespace Amanzi
{

Transient_Richards_PK::Transient_Richards_PK(Teuchos::ParameterList &plist, const Teuchos::RCP<const Flow_State> FS_) : FS(FS_), richards_plist(plist)
{
  // Add some parameters to the Richards problem constructor parameter list.
  Teuchos::ParameterList &rp_list = plist.sublist("Richards Problem");
  rp_list.set("fluid density", FS->fluid_density());
  rp_list.set("fluid viscosity", FS->fluid_viscosity());
  const double *gravity = FS->gravity();
  //TODO: assuming gravity[0] = gravity[1] = 0 -- needs to be reconciled somehow
  rp_list.set("gravity", -gravity[2]);
  
  // Create the Richards flow problem.
  Teuchos::ParameterList rlist = richards_plist.sublist("Richards Problem");
  problem = new RichardsProblem(FS->mesh(), rlist);

  // ss_t0 = rlist.get<double>("Steady state calculation initial time");
  // ss_t1 = rlist.get<double>("Steady state calculation final time");
  // ss_h0 = rlist.get<double>("Steady state calculation initial time step");
  // ss_z =  rlist.get<double>("Steady state calculation initial hydrostatic pressure height");

  // Create the solution vectors.
  solution = new Epetra_Vector(problem->Map());
  pressure_cells = problem->CreateCellView(*solution);
  pressure_faces = problem->CreateFaceView(*solution);
  richards_flux = new Epetra_Vector(problem->FaceMap());

  // create the time stepper...

  // first the Richards model evaluator
  Teuchos::ParameterList &rme_list = rlist.sublist("Richards model evaluator");
  RME = new RichardsModelEvaluator(problem, rme_list, problem->Map(), FS);  

  // then the BDF2 solver
  Teuchos::RCP<Teuchos::ParameterList> bdf2_list_p(new Teuchos::ParameterList(rlist.sublist("Time integrator")));

  time_stepper = new BDF2::Dae(*RME, problem->Map());
  time_stepper->setParameterList(bdf2_list_p);

};

Transient_Richards_PK::~Transient_Richards_PK()
{
  delete richards_flux;
  delete pressure_cells;
  delete pressure_faces;
  delete solution;
  delete problem;
};


int Transient_Richards_PK::advance_to_steady_state()
{
  // Set problem parameters.
  problem->SetPermeability(FS->permeability());
  problem->SetFlowState(FS);

  double t0 = ss_t0;
  double t1 = ss_t1;
  double h =  ss_h0;
  double hnext;

  // create udot

  problem->SetInitialPressureProfileCells(ss_z,pressure_cells);
  problem->SetInitialPressureProfileFaces(ss_z,pressure_faces);

  Epetra_Vector udot(problem->Map());
  problem->Compute_udot(t0,  *solution, udot);

  time_stepper->set_initial_state(t0, *solution, udot);

  int errc;
  RME->update_precon(t0, *solution, h, errc);

  // iterate
  int i = 0;
  double tlast = t0;

  do {
    
    time_stepper->bdf2_step(h,0.0,*solution,hnext);

    time_stepper->commit_solution(h,*solution);

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

int Transient_Richards_PK::init_transient(double t0, double h_)
{
  h = h_;
  hnext = h_;

  // Set problem parameters.
  problem->SetPermeability(FS->permeability());
  problem->SetFlowState(FS);

  Epetra_Vector udot(problem->Map());
  problem->Compute_udot(t0,  *solution, udot);
  
  time_stepper->set_initial_state(t0, *solution, udot);
  
  int errc;
  RME->update_precon(t0, *solution, h, errc);
}


int Transient_Richards_PK::advance_transient(double h) 
{
  // Set problem parameters.
  problem->SetPermeability(FS->permeability());
  problem->SetFlowState(FS);

  time_stepper->bdf2_step(h,0.0,*solution,hnext);
  time_stepper->commit_solution(h,*solution);  

  time_stepper->write_bdf2_stepping_statistics();
}






void Transient_Richards_PK::GetSaturation(Epetra_Vector &s) const
{
  //for (int i = 0; i < s.MyLength(); ++i) s[i] = 1.0;

  problem->DeriveVanGenuchtenSaturation(*pressure_cells, s);

}


} // close namespace Amanzi
