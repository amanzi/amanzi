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
#include "State.hpp"

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
    mesh = mesh_factory(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 4, 4, 4);
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


  // create the state
  Teuchos::ParameterList state_plist;
  state_plist.set<int>("Number of mesh blocks",1);
  state_plist.set<int>("Number of component concentrations",1);
  state_plist.set<double>("Constant water density",1000.0);
  state_plist.set<double>("Constant water saturation",1.0);
  state_plist.set<double>("Constant viscosity",1.0);
  state_plist.set<double>("Gravity x",0.0);
  state_plist.set<double>("Gravity y",0.0);
  state_plist.set<double>("Gravity z",0.0);
  
  Teuchos::ParameterList& mb1_plist = state_plist.sublist("Mesh block 1");
  mb1_plist.set<int>("Mesh block ID",0);
  mb1_plist.set<double>("Constant porosity",0.2);
  mb1_plist.set<double>("Constant permeability",10.0);
  mb1_plist.set<double>("Constant component concentration 0",0.0);
  mb1_plist.set<double>("Constant Darcy flux x",0.0);
  mb1_plist.set<double>("Constant Darcy flux y",0.0);
  mb1_plist.set<double>("Constant Darcy flux z",0.0);


  Teuchos::RCP<State> S = Teuchos::rcp(new State(state_plist, mesh));
  Teuchos::RCP<Flow_State> FS = Teuchos::rcp(new Flow_State(S));

  
  // PROBLEM
  Teuchos::ParameterList pl;
  pl.set("van Genuchten m",0.5);
  pl.set("van Genuchten alpha",0.00005);
  pl.set("van Genuchten residual saturation",0.40969);
  pl.set("atmospheric pressure", 0.0);
  RichardsProblem problem(mesh, pl, bc);

  // MODEL PARAMETERS
  problem.SetFluidDensity(FS->fluid_density());
  problem.SetFluidViscosity(FS->fluid_viscosity());
  problem.SetPermeability(1.0);
  problem.SetGravity(FS->gravity());

  // Create the BDF2 parameter list.
  Teuchos::RCP<Teuchos::ParameterList> bdf2_param_p(new Teuchos::ParameterList);
  Teuchos::ParameterList &bdf2_param = *bdf2_param_p.get();

  // BDF2 Paramters
  bdf2_param.set("Nonlinear solver max iterations",20);
  bdf2_param.set("Nonlinear solver tolerance",0.01);
  bdf2_param.set("NKA max vectors", 10);
  bdf2_param.set("NKA drop tolerance", 5e-2);

  // set the BDF2 verbosity level
  bdf2_param.sublist("VerboseObject").set("Verbosity Level","high");



  // create the Richards Model Evaluator
  Teuchos::ParameterList plist;
  plist.sublist("VerboseObject").set("Verbosity Level","none");
  
  
  RichardsModelEvaluator RME(&problem, plist, problem.Map(), FS);

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
 
  
  

  
  // Write a GMV file with the cell pressures.
  std::string gmv_file("pressure.gmv");
  GMV::open_data_file(*mesh, gmv_file);
  GMV::start_data();
  Epetra_Vector &Pcell = *(problem.CreateCellView(u));
  GMV::write_cell_data(Pcell, std::string("pressure"));
  GMV::close_data_file();
  delete &Pcell;
  } while (t1 >= tlast);

  return 0;
}
