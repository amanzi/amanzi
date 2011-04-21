#include "Mesh_maps_moab.hh"

#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Epetra_BlockMap.h"
#include "Epetra_MpiComm.h"
#include "mpi.h"

#include "MeshFactory.hh"
#include "FrameworkTraits.hh"

#include "RichardsProblem.hpp"
#include "RichardsModelEvaluator.hpp"
#include "BDF2_Dae.hpp"
#include "gmv_mesh.hh"

int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  Epetra_MpiComm comm(MPI_COMM_WORLD);

  // MESH
  Mesh::MeshFactory mesh_factory(comm);
  Teuchos::RCP<Mesh_maps_base> mesh;
  Mesh::FrameworkPreference pref;
  
  if (framework_available(Mesh::Simple)) {
    pref.clear(); pref.push_back(Mesh::Simple);
    mesh_factory.preference(pref);
    mesh = mesh_factory(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 3, 3, 3);
  }

  // BOUNDARY CONDITIONS
  Teuchos::ParameterList bc_params;
  bc_params.set("number of BCs", 6);
  Teuchos::ParameterList &bc_left = bc_params.sublist("BC00");
  bc_left.set("Side set ID", 1);
  bc_left.set("Type", "Pressure Constant");
  bc_left.set("BC value", 1.0);
  Teuchos::ParameterList &bc_right = bc_params.sublist("BC01");
  bc_right.set("Side set ID", 3);
  bc_right.set("Type", "Pressure Constant");
  bc_right.set("BC value", 0.0);
  Teuchos::ParameterList &bc_front = bc_params.sublist("BC02");
  bc_front.set("Side set ID", 0);
  bc_front.set("Type", "No Flow");
  Teuchos::ParameterList &bc_back = bc_params.sublist("BC03");
  bc_back.set("Side set ID", 2);
  bc_back.set("Type", "No Flow");
  Teuchos::ParameterList &bc_bottom = bc_params.sublist("BC04");
  bc_bottom.set("Side set ID", 4);
  bc_bottom.set("Type", "No Flow");
  Teuchos::ParameterList &bc_top = bc_params.sublist("BC05");
  bc_top.set("Side set ID", 5);
  bc_top.set("Type", "No Flow");
  Teuchos::RCP<FlowBC> bc(new FlowBC(bc_params, mesh));
  
  // PROBLEM
  Teuchos::ParameterList pl;
  pl.set("van Genuchten m",0.5);
  pl.set("van Genuchten alpha",0.00005);
  pl.set("van Genuchten residual saturation",0.40969);
  pl.set("atmospheric pressure", 0.0);
  RichardsProblem problem(mesh, pl, bc);

  // MODEL PARAMETERS
  problem.SetFluidDensity(1.0);
  problem.SetFluidViscosity(1.0);
  problem.SetPermeability(1.0);
  double g[3] = { 0.0 };
  g[2] = -9.8;
  problem.SetGravity(g);

  // Create the BDF2 parameter list.
  Teuchos::RCP<Teuchos::ParameterList> bdf2_param_p(new Teuchos::ParameterList);
  Teuchos::ParameterList &bdf2_param = *bdf2_param_p.get();

  // BDF2 Paramters
  bdf2_param.set("Nonlinear solver max iterations",10);
  bdf2_param.set("Nonlinear solver tolerance",0.01);
  bdf2_param.set("NKA max vectors", 4);
  bdf2_param.set("NKA drop tolerance", 5e-2);

  // set the BDF2 verbosity level
  bdf2_param.sublist("VerboseObject").set("Verbosity Level","high");

  // create the Richards Model Evaluator
  Teuchos::ParameterList plist;
  plist.sublist("VerboseObject").set("Verbosity Level","none");
  
  RichardsModelEvaluator RME(&problem, plist, problem.Map());

  // create the time stepping object
  BDF2::Dae TS( RME, problem.Map() );
  TS.setParameterList(bdf2_param_p);

  // create the initial condition
  Epetra_Vector u(problem.Map());
  u.PutScalar(0.0);

  // set intial and final time
  double t0 = 0.0;
  double t1 = 1.0;

  // compute the initial udot
  Epetra_Vector udot(problem.Map());
  problem.Compute_udot(t0,u,udot);
  
  // initial time step
  double h = 1.0e-5;
  double hnext;

  // intialize the state of the time stepper
  TS.set_initial_state(t0,u,udot);
 
  int errc;
  RME.update_precon(t0, u, h, errc);

  // iterate
  int i = 0;
  double tlast;
  do {
    
    TS.bdf2_step(h,0.0,20,u,hnext);
    
    TS.commit_solution(h,u);
    
    TS.write_bdf2_stepping_statistics();

    h = hnext;
    i++;

    tlast=TS.most_recent_time();
  } while (t1 >= tlast);
 
  
  

  
  // Write a GMV file with the cell pressures.
  std::string gmv_file("pressure.gmv");
  GMV::open_data_file(*mesh, gmv_file);
  GMV::start_data();
  Epetra_Vector &Pcell = *(problem.CreateCellView(u));
  GMV::write_cell_data(Pcell, std::string("pressure"));
  GMV::close_data_file();
  delete &Pcell;

  return 0;
}
