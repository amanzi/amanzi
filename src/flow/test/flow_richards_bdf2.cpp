/*
The flow component of the Amanzi code, richards unit tests.
License: BSD
Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>

#include "UnitTest++.h"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "Mesh_simple.hh"
#include "State.hpp"
#include "Flow_State.hpp"
#include "RichardsProblem.hpp"
#include "RichardsModelEvaluator.hpp"
#include "BDF2_Dae.hpp"

#include "gmv_mesh.hh"
#include "cgns_mesh_par.hh"
#include "cgns_mesh.hh"


/* **************************************************************** */
TEST(FLOW_RICHARDS) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;

cout << "Test: 3D Richards on a cubic mesh" << endl;
#ifdef HAVE_MPI
  Epetra_MpiComm  *comm = new Epetra_MpiComm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm  *comm = new Epetra_SerialComm();
#endif

  /* read parameter list */
  ParameterList parameter_list;
  string xmlFileName = "test/flow_richards_bdf2.xml";
  updateParametersFromXmlFile(xmlFileName, &parameter_list);

  // create an SIMPLE mesh framework 
  ParameterList region_list = parameter_list.get<Teuchos::ParameterList>("Regions");
  GeometricModelPtr gm = new GeometricModel(3, region_list);
  RCP<Mesh> mesh = rcp(new Mesh_simple(0.0,0.0,-68.0, 1.0,1.0,0.0, 1, 1, 136, comm, gm)); 

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
  Teuchos::RCP<Amanzi::Flow_State> FS = Teuchos::rcp(new Amanzi::Flow_State(S));

  
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

  Amanzi::RichardsProblem problem(mesh, vGl, bc);

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
  
  
  Amanzi::RichardsModelEvaluator RME(&problem, plist, problem.Map(), FS);

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

  Amanzi::CGNS::create_mesh_file(*mesh, cgns_filename);

  Amanzi::CGNS::open_data_file(cgns_filename);
  Amanzi::CGNS::create_timestep(0.0, 0, Amanzi::AmanziMesh::CELL);

  Amanzi::CGNS::write_field_data(*(S->get_pressure()), "pressure");
  Amanzi::CGNS::write_field_data(*(S->get_permeability()), "permeability");
  
  // iterate
  int i = 0;
  double tlast = t0;
  do {    
    TS.bdf2_step(h,0.0,u,hnext);
    
    TS.commit_solution(h,u);

    S->advance_time(h);

    // update the state, but only the cell values of pressure
    S->update_pressure( * problem.CreateCellView(u) );

    TS.write_bdf2_stepping_statistics();

    h = hnext;
    i++;

    tlast=TS.most_recent_time();

    if ( i%5 == 1 ) { 
      Amanzi::CGNS::open_data_file(cgns_filename);
      Amanzi::CGNS::create_timestep(S->get_time(), i, Amanzi::AmanziMesh::CELL);

      Amanzi::CGNS::write_field_data(*(S->get_pressure()), "pressure");
      Amanzi::CGNS::write_field_data(*(S->get_permeability()), "permeability");
    }
  
  } while (t1 >= tlast);

  Amanzi::CGNS::open_data_file(cgns_filename);
  Amanzi::CGNS::create_timestep(S->get_time(), i, Amanzi::AmanziMesh::CELL);
  
  Amanzi::CGNS::write_field_data(*(S->get_pressure()), "pressure");
  Amanzi::CGNS::write_field_data(*(S->get_permeability()), "permeability");
 
  delete comm;
}
