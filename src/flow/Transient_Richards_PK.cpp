#include "Teuchos_VerboseObjectParameterListHelpers.hpp"

#include "Mesh.hh"
#include "Transient_Richards_PK.hpp"
#include "RichardsProblem.hpp"

namespace Amanzi {

Transient_Richards_PK::Transient_Richards_PK(Teuchos::ParameterList &plist_, 
                                             const Teuchos::RCP<Flow_State> FS_) : FS(FS_), plist(plist_)
{
  // set the line prefix for output
  this->setLinePrefix("Amanzi::Transient_Richards_PK ");
  // make sure that the line prefix is printed
  this->getOStream()->setShowLinePrefix(true);

  // Read the sublist for verbosity settings.
  Teuchos::readVerboseObjectSublist(&plist,this);

  using Teuchos::OSTab;
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  OSTab tab = this->getOSTab(); // This sets the line prefix and adds one tab

  // Add some parameters to the Richards problem constructor parameter list.
  Teuchos::ParameterList &rp_list = plist.sublist("Richards Problem");
  rp_list.set("fluid density", FS->get_fluid_density());
  rp_list.set("fluid viscosity", FS->get_fluid_viscosity());
  const double *gravity = FS->get_gravity();
  //TODO: assuming gravity[0] = gravity[1] = 0 -- needs to be reconciled somehow
  rp_list.set("gravity", -gravity[2]);
  
  // Create the Richards flow problem.
  problem = new RichardsProblem(FS->get_mesh_maps(), rp_list);

  // Create the solution vectors.
  solution = new Epetra_Vector(problem->Map());
  pressure_cells = problem->CreateCellView(*solution);
  pressure_faces = problem->CreateFaceView(*solution);
  richards_flux = new Epetra_Vector(problem->FaceMap());
  
  // get the pressure from the flow state
  *pressure_cells = FS->get_pressure();
  // and compute approximate face pressures
  approximate_face_pressure(*pressure_cells, *pressure_faces);

  // first the Richards model evaluator
  Teuchos::ParameterList &rme_list = rp_list.sublist("Richards model evaluator");
  RME = new RichardsModelEvaluator(problem, rme_list, problem->Map(), FS);  

  // then the BDF2 solver
  Teuchos::RCP<Teuchos::ParameterList> bdf2_list_p(new Teuchos::ParameterList(rp_list.sublist("Time integrator")));

  time_stepper = new BDF2::Dae(*RME, problem->Map());
  time_stepper->setParameterList(bdf2_list_p);


  // then the BDF1 solver
  Teuchos::RCP<Teuchos::ParameterList> bdf1_list_p(new Teuchos::ParameterList(rp_list.sublist("steady time integrator")));

  steady_time_stepper = new Amanzi::BDF1Dae(*RME, problem->Map());

  steady_time_stepper->setParameterList(bdf1_list_p);
  
  // initialize the water saturation for vis
  GetSaturation( FS->get_water_saturation() );

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
  problem->set_absolute_permeability(FS->get_vertical_permeability(), FS->get_horizontal_permeability());
  problem->set_flow_state(FS);

  double t0 = ss_t0;
  double t1 = ss_t1;
  double h =  ss_h0;
  double hnext;

  // create udot

  problem->set_pressure_cells(ss_z, pressure_cells);
  problem->set_pressure_faces(ss_z, pressure_faces);

  Epetra_Vector udot(problem->Map());
  problem->compute_udot(t0, *solution, udot);

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
  problem->set_absolute_permeability(FS->get_vertical_permeability(), FS->get_horizontal_permeability());
  problem->set_flow_state(FS);

  Epetra_Vector udot(problem->Map());
  problem->compute_udot(t0, *solution, udot);
  
  time_stepper->set_initial_state(t0, *solution, udot);
  
  int errc;
  RME->update_precon(t0, *solution, h, errc);
  RME->update_norm(0.001,1.0);
}

int Transient_Richards_PK::init_steady(double t0, double h_)
{
  h = h_;
  hnext = h_;

  // Set problem parameters.
  problem->set_absolute_permeability(FS->get_vertical_permeability(), FS->get_horizontal_permeability());
  problem->set_flow_state(FS);

  Epetra_Vector udot(problem->Map());
  problem->compute_udot(t0, *solution, udot);
  
  steady_time_stepper->set_initial_state(t0, *solution, udot);
  
  int errc;
  RME->update_precon(t0, *solution, h, errc);
  RME->update_norm(0.0, 1.0); // we run the steady calculation with just an absolute norm
}

int Transient_Richards_PK::advance_transient(double h) 
{
  //using Teuchos::OSTab;
  //Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  //Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  //OSTab tab = this->getOSTab(); // This sets the line prefix and adds one tab

  // Set problem parameters
  problem->set_absolute_permeability(FS->get_vertical_permeability(), FS->get_horizontal_permeability());
  problem->set_flow_state(FS);

  time_stepper->bdf2_step(h,0.0,*solution,hnext);
  time_stepper->commit_solution(h,*solution);  

  time_stepper->write_bdf2_stepping_statistics();

  
  //FS->get_pressure() = *pressure_cells;
  //FS->get_prev_water_saturation() = FS->get_water_saturation();
  //GetSaturation( FS->get_water_saturation() );

  //double l1_error;
  //problem->DeriveDarcyFlux(*solution, *richards_flux, l1_error);
  //if (out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_HIGH) ) {
  //  *out << "L1 norm of the Richards flux discrepancy = " << l1_error << std::endl;
 // }
}


int Transient_Richards_PK::advance_steady(double h) 
{
  //using Teuchos::OSTab;
  //Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  //Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  //OSTab tab = this->getOSTab(); // This sets the line prefix and adds one tab

  // Set problem parameters
  problem->set_absolute_permeability(FS->get_vertical_permeability(), FS->get_horizontal_permeability());
  problem->set_flow_state(FS);

  steady_time_stepper->bdf1_step(h,*solution,hnext);
  steady_time_stepper->commit_solution(h,*solution);  

  steady_time_stepper->write_bdf1_stepping_statistics();

  //FS->get_pressure() = *pressure_cells;
  //FS->get_prev_water_saturation() = FS->get_water_saturation();
  //GetSaturation( FS->get_water_saturation() );  

  //double l1_error;
  //problem->DeriveDarcyFlux(*solution, *richards_flux, l1_error);
  //if (out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_HIGH) ) {
  //  *out << "L1 norm of the Richards flux discrepancy = " << l1_error << std::endl;
 // }
}




void Transient_Richards_PK::GetSaturation(Epetra_Vector &s) const
{
  problem->DeriveVanGenuchtenSaturation(*pressure_cells, s);
}
  
void  Transient_Richards_PK::commit_state(Teuchos::RCP<Flow_State> FS) 
{
  using Teuchos::OSTab;
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  OSTab tab = this->getOSTab(); // This sets the line prefix and adds one tab

  FS->get_pressure() = *pressure_cells;
  FS->get_prev_water_saturation() = FS->get_water_saturation();
  GetSaturation( FS->get_water_saturation() );  

  double l1_error;
  problem->DeriveDarcyFlux(*solution, *richards_flux, l1_error);
  FS->get_darcy_flux() = *richards_flux;
  if (out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_HIGH) ) {
    *out << "L1 norm of the Richards flux discrepancy = " << l1_error << std::endl;
  }

  // derive the velocity vector
  GetVelocity(FS->get_darcy_velocity());


  // now compute the discrepancy between the current and previous saturations
  if (out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_HIGH) ) {
    Epetra_Vector s_diff (FS->get_water_saturation());
    s_diff = FS->get_water_saturation();
    s_diff.Update(-1.0,FS->get_prev_water_saturation(),1.0);
    double norm;
    s_diff.NormInf(&norm);
    *out << "discrepancy between current and previous saturation = " << norm << std::endl;
  }
}


void Transient_Richards_PK::approximate_face_pressure(const Epetra_Vector& cell_pressure, Epetra_Vector& face_pressure)
{
  // make a vector that includes the overlap region
  Epetra_Vector cell_pressure_ovl(FS->get_mesh_maps()->cell_epetra_map(true) );

  // import the cell pressures to this vector to get access to the overlap
  Epetra_Import cell_importer(FS->get_mesh_maps()->cell_epetra_map(true),FS->get_mesh_maps()->cell_epetra_map(false));
  
  cell_pressure_ovl.Import(cell_pressure, cell_importer, Insert);

  // loop over all faces
  for (int iface=0; iface<face_pressure.MyLength(); iface++)
    {
      // find the neighbor cell of the current face
      Amanzi::AmanziMesh::Entity_ID_List cells;
      FS->get_mesh_maps()->face_get_cells(iface,Amanzi::AmanziMesh::USED, &cells);
      
      double cp(0.0);
      for (unsigned int it=0; it<cells.size(); it++)
	{
	  cp += cell_pressure_ovl[cells[it]];
	}
      // compute an average pressure for the face
      face_pressure[iface] = cp / cells.size();
    }
}



}  // close namespace Amanzi
