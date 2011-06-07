#include "Mesh_maps_simple.hh"

#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Epetra_BlockMap.h"
#include "Epetra_MpiComm.h"
#include "mpi.h"

//#include "MeshFactory.hh"
//#include "FrameworkTraits.hh"


#include "RichardsProblem.hpp"
#include "RichardsModelEvaluator.hpp"
#include "BDF2_Dae.hpp"
#include "gmv_mesh.hh"
#include "State.hpp"

#include "cgns_mesh_par.hh"
#include "cgns_mesh.hh"

int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  Epetra_MpiComm comm(MPI_COMM_WORLD);

  // simple mesh with two layers
  Teuchos::ParameterList smeshlist;
  smeshlist.set<int>("Numer of Cells in X",1);
  smeshlist.set<int>("Numer of Cells in Y",1);
  smeshlist.set<int>("Numer of Cells in Z",136);
  smeshlist.set<double>("X_Min",0.0);
  smeshlist.set<double>("X_Max",1.0);
  smeshlist.set<double>("Y_Min",0.0);
  smeshlist.set<double>("Y_Max",1.0);
  smeshlist.set<double>("Z_Min",-68.0);
  smeshlist.set<double>("Z_Max",0.0);

  smeshlist.set<int>("Number of mesh blocks",5);

  {
    Teuchos::ParameterList &mblist = smeshlist.sublist("Mesh block 1");
    mblist.set<double>("Z0",-69.0);
    mblist.set<double>("Z1",-48.0);
  }
  {
    Teuchos::ParameterList &mblist = smeshlist.sublist("Mesh block 2");
    mblist.set<double>("Z0",-48.0);
    mblist.set<double>("Z1",-45.0);
  }
  {
    Teuchos::ParameterList &mblist = smeshlist.sublist("Mesh block 3");
    mblist.set<double>("Z0",-45.0);
    mblist.set<double>("Z1",-39.0);
  }  
  {
    Teuchos::ParameterList &mblist = smeshlist.sublist("Mesh block 4");
    mblist.set<double>("Z0",-39.0);
    mblist.set<double>("Z1",-16.0);
  }
  {
    Teuchos::ParameterList &mblist = smeshlist.sublist("Mesh block 5");
    mblist.set<double>("Z0",-16.0);
    mblist.set<double>("Z1",1.0);
  }

  Teuchos::RCP<Mesh_maps_base> mesh = Teuchos::rcp(new Mesh_maps_simple(smeshlist,&comm));

  // BOUNDARY CONDITIONS
  Teuchos::ParameterList bc_params;
  bc_params.set("number of BCs", 6);
  Teuchos::ParameterList &bc_left = bc_params.sublist("BC00");
  bc_left.set("Side set ID", 1);
  bc_left.set("Type", "No Flow");
  Teuchos::ParameterList &bc_right = bc_params.sublist("BC01");
  bc_right.set("Side set ID", 3);
  bc_right.set("Type", "No Flow");
  Teuchos::ParameterList &bc_front = bc_params.sublist("BC02");
  bc_front.set("Side set ID", 0);
  bc_front.set("Type", "No Flow");
  Teuchos::ParameterList &bc_back = bc_params.sublist("BC03");
  bc_back.set("Side set ID", 2);
  bc_back.set("Type", "No Flow");
  Teuchos::ParameterList &bc_bottom = bc_params.sublist("BC04");
  bc_bottom.set("Side set ID", 4);
  bc_bottom.set("Type", "Pressure Constant");
  bc_bottom.set("BC value", 160125.0);  
  Teuchos::ParameterList &bc_top = bc_params.sublist("BC05");
  bc_top.set("Side set ID", 5);
  bc_top.set("Type", "Pressure Constant");
  bc_top.set("BC value", -506275.0);   
 
  Teuchos::RCP<FlowBC> bc(new FlowBC(bc_params, mesh));


  // create the state
  Teuchos::ParameterList state_plist;
  state_plist.set<int>("Number of component concentrations",1);
  state_plist.set<double>("Constant water density",1000.0);
  state_plist.set<double>("Constant water saturation",1.0);
  state_plist.set<double>("Constant viscosity",1.0);
  state_plist.set<double>("Gravity x",0.0);
  state_plist.set<double>("Gravity y",0.0);
  state_plist.set<double>("Gravity z",-9.8);

  state_plist.set<int>("Number of mesh blocks",5);

  {
    Teuchos::ParameterList& mb_plist = state_plist.sublist("Mesh block 1");
    mb_plist.set<int>("Mesh block ID",0);
    mb_plist.set<double>("Constant porosity",0.1643);
    mb_plist.set<double>("Constant permeability",7.33e-14);
    mb_plist.set<double>("Constant component concentration 0",0.0);
    mb_plist.set<double>("Constant Darcy flux x",0.0);
    mb_plist.set<double>("Constant Darcy flux y",0.0);
    mb_plist.set<double>("Constant Darcy flux z",0.0);
  }
  {
    Teuchos::ParameterList& mb_plist = state_plist.sublist("Mesh block 2");
    mb_plist.set<int>("Mesh block ID",1);
    mb_plist.set<double>("Constant porosity",0.2625);
    mb_plist.set<double>("Constant permeability",5.27e-14);
    mb_plist.set<double>("Constant component concentration 0",0.0);
    mb_plist.set<double>("Constant Darcy flux x",0.0);
    mb_plist.set<double>("Constant Darcy flux y",0.0);
    mb_plist.set<double>("Constant Darcy flux z",0.0);
  }    
  {
    Teuchos::ParameterList& mb_plist = state_plist.sublist("Mesh block 3");
    mb_plist.set<int>("Mesh block ID",2);
    mb_plist.set<double>("Constant porosity",0.4223);
    mb_plist.set<double>("Constant permeability",1.37e-14);
    mb_plist.set<double>("Constant component concentration 0",0.0);
    mb_plist.set<double>("Constant Darcy flux x",0.0);
    mb_plist.set<double>("Constant Darcy flux y",0.0);
    mb_plist.set<double>("Constant Darcy flux z",0.0);
  }
  {
    Teuchos::ParameterList& mb_plist = state_plist.sublist("Mesh block 4");
    mb_plist.set<int>("Mesh block ID",3);
    mb_plist.set<double>("Constant porosity",0.3586);
    mb_plist.set<double>("Constant permeability",1.23e-13);
    mb_plist.set<double>("Constant component concentration 0",0.0);
    mb_plist.set<double>("Constant Darcy flux x",0.0);
    mb_plist.set<double>("Constant Darcy flux y",0.0);
    mb_plist.set<double>("Constant Darcy flux z",0.0);
  }    
  {
    Teuchos::ParameterList& mb_plist = state_plist.sublist("Mesh block 5");
    mb_plist.set<int>("Mesh block ID",4);
    mb_plist.set<double>("Constant porosity",0.2585);
    mb_plist.set<double>("Constant permeability",1.24e-12);
    mb_plist.set<double>("Constant component concentration 0",0.0);
    mb_plist.set<double>("Constant Darcy flux x",0.0);
    mb_plist.set<double>("Constant Darcy flux y",0.0);
    mb_plist.set<double>("Constant Darcy flux z",0.0);
  }    

  Teuchos::RCP<State> S = Teuchos::rcp(new State(state_plist, mesh));
  Teuchos::RCP<Flow_State> FS = Teuchos::rcp(new Flow_State(S));

  
  // Water Retention Model (here we use van Genuchten)
  Teuchos::ParameterList vGl;
  vGl.set<double>("Atmospheric pressure", 101325.0);

  Teuchos::ParameterList& vGsl = vGl.sublist("Water retention models");

  {
    Teuchos::ParameterList& mb_vGl = vGsl.sublist("WRM 0");
    mb_vGl.set<string>("Water retention model","van Genuchten");
    mb_vGl.set<int>("Region ID", 0);
    mb_vGl.set<double>("van Genuchten m",0.392);
    mb_vGl.set<double>("van Genuchten alpha",6.33e-5);
    mb_vGl.set<double>("van Genuchten residual saturation",0.0609);
  }
  {
    Teuchos::ParameterList& mb_vGl = vGsl.sublist("WRM 1");
    mb_vGl.set<string>("Water retention model","van Genuchten");
    mb_vGl.set<int>("Region ID", 1);
    mb_vGl.set<double>("van Genuchten m",0.386);
    mb_vGl.set<double>("van Genuchten alpha",2.961e-5);
    mb_vGl.set<double>("van Genuchten residual saturation",0.213);
  }    
  {
    Teuchos::ParameterList& mb_vGl = vGsl.sublist("WRM 2");
    mb_vGl.set<string>("Water retention model","van Genuchten");
    mb_vGl.set<int>("Region ID", 2);
    mb_vGl.set<double>("van Genuchten m",0.456);
    mb_vGl.set<double>("van Genuchten alpha",6.84e-5);
    mb_vGl.set<double>("van Genuchten residual saturation",0.2595);
  }
  {
    Teuchos::ParameterList& mb_vGl = vGsl.sublist("WRM 3");
    mb_vGl.set<string>("Water retention model","van Genuchten");
    mb_vGl.set<int>("Region ID", 3);
    mb_vGl.set<double>("van Genuchten m",0.469);
    mb_vGl.set<double>("van Genuchten alpha",9.39e-5);
    mb_vGl.set<double>("van Genuchten residual saturation",0.0837);
  }    
  {
    Teuchos::ParameterList& mb_vGl = vGsl.sublist("WRM 4");
    mb_vGl.set<string>("Water retention model","van Genuchten");
    mb_vGl.set<int>("Region ID", 4);
    mb_vGl.set<double>("van Genuchten m",0.658);
    mb_vGl.set<double>("van Genuchten alpha",1.008e-3);
    mb_vGl.set<double>("van Genuchten residual saturation",0.0774);
  }    

  RichardsProblem problem(mesh, vGl, bc);

  // MODEL PARAMETERS
  problem.SetFluidDensity(FS->fluid_density());
  problem.SetFluidViscosity(FS->fluid_viscosity());
  problem.SetPermeability(1.0);
  problem.SetGravity(FS->gravity());
  problem.SetFlowState(FS);

  // Create the BDF2 parameter list.
  Teuchos::RCP<Teuchos::ParameterList> bdf2_param_p(new Teuchos::ParameterList);
  Teuchos::ParameterList &bdf2_param = *bdf2_param_p.get();

  // BDF2 Paramters
  bdf2_param.set("Nonlinear solver max iterations",10);
  bdf2_param.set("Nonlinear solver tolerance",0.01);
  bdf2_param.set("NKA max vectors", 5);
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

  // first for the cells
  Epetra_Vector * ucells =  problem.CreateCellView(u);
  for (int j=0; j<ucells->MyLength(); j++)
    {
      std::vector<double> coords;
      coords.resize(24);
      mesh->cell_to_coordinates(j, coords.begin(), coords.end());

      // average the x coordinates
      double zavg = 0.0;
      for (int k=2; k<24; k+=3)
	zavg += coords[k];
      zavg /= 8.0;

      (*ucells)[j] = 101325.0 - 9800 * (zavg + 62.0);
    }

  // then for faces
  Epetra_Vector * ufaces =  problem.CreateFaceView(u);
  for (int j=0; j<ufaces->MyLength(); j++)
    {
      std::vector<double> coords;
      coords.resize(12);
      mesh->face_to_coordinates(j, coords.begin(), coords.end());

      // average the x coordinates
      double zavg = 0.0;
      for (int k=2; k<12; k+=3)
	zavg += coords[k];
      zavg /= 4.0;

      (*ufaces)[j] = 101325.0 - 9800 * (zavg + 62.0);
    }
  
  // initialize the state
  S->update_pressure( * problem.CreateCellView(u) );
  S->set_time(0.0);

  // set intial and final time
  double t0 = 0.0;
  double t1 = 1000.0;

  // compute the initial udot
  Epetra_Vector udot(problem.Map());
  problem.Compute_udot(t0,u,udot);

  udot.PutScalar(0.0);

  // initial time step
  double h = 1.0e-5;
  double hnext;

  // intialize the state of the time stepper
  TS.set_initial_state(t0,u,udot);
 
  int errc;
  RME.update_precon(t0, u, h, errc);


  // set up output
  std::string cgns_filename = "out.cgns";

  CGNS::create_mesh_file(*mesh, cgns_filename);

  CGNS::open_data_file(cgns_filename);
  CGNS::create_timestep(0.0, 0, Mesh_data::CELL);

  CGNS::write_field_data(*(S->get_pressure()), "pressure");
  CGNS::write_field_data(*(S->get_permeability()), "permeability");
  
  // iterate
  int i = 0;
  double tlast = t0;
  do {
    
    TS.bdf2_step(h,0.0,20,u,hnext);
    
    TS.commit_solution(h,u);

    S->advance_time(h);

    // update the state, but only the cell values of pressure
    S->update_pressure( * problem.CreateCellView(u) );

    TS.write_bdf2_stepping_statistics();

    h = hnext;
    i++;

    tlast=TS.most_recent_time();
    

    if ( i%5 == 1 ) 
      { 
	CGNS::open_data_file(cgns_filename);
	CGNS::create_timestep(S->get_time(), i, Mesh_data::CELL);

	CGNS::write_field_data(*(S->get_pressure()), "pressure");
	CGNS::write_field_data(*(S->get_permeability()), "permeability");

      }

  
  } while (t1 >= tlast);

  CGNS::open_data_file(cgns_filename);
  CGNS::create_timestep(S->get_time(), i, Mesh_data::CELL);
  
  CGNS::write_field_data(*(S->get_pressure()), "pressure");
  CGNS::write_field_data(*(S->get_permeability()), "permeability");

  

  return 0;
}
