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
TEST(FLOW_1D_RICHARDS) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;

cout << "Test: 1D Richards, 2-layer model" << endl;
#ifdef HAVE_MPI
  Epetra_MpiComm  *comm = new Epetra_MpiComm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm  *comm = new Epetra_SerialComm();
#endif

  /* read parameter list */
  ParameterList parameter_list;
  string xmlFileName = "test/flow_richards_bdf2.xml"; //new-bc-cribs-1D.xml";
  updateParametersFromXmlFile(xmlFileName, &parameter_list);

  // create an SIMPLE mesh framework 
  ParameterList region_list = parameter_list.get<Teuchos::ParameterList>("Regions");
  GeometricModelPtr gm = new GeometricModel(3, region_list, (Epetra_MpiComm *)comm);
  RCP<Mesh> mesh = rcp(new Mesh_simple(0.0,0.0,-10.0, 1.0,1.0,0.0, 1, 1, 100, comm, gm)); 

  // create the state
  ParameterList state_list = parameter_list.get<Teuchos::ParameterList>("State");
  RCP<State> S = Teuchos::rcp(new State(state_list, mesh));
  RCP<Amanzi::Flow_State> FS = Teuchos::rcp(new Amanzi::Flow_State(S));

  // create Richards problem
  ParameterList flow_list = parameter_list.get<Teuchos::ParameterList>("Flow");
  ParameterList richards_list = flow_list.get<Teuchos::ParameterList>("Richards Problem");

  richards_list.set("fluid density", FS->get_fluid_density());  // will be removed in new version
  richards_list.set("fluid viscosity", FS->get_fluid_viscosity());
  const double *gravity = FS->get_gravity();
  richards_list.set("gravity", gravity[2]);

  RichardsProblem problem(mesh, richards_list);
  problem.set_flow_state(FS);
  problem.set_absolute_permeability(FS->get_vertical_permeability());

  // create the Richards Model Evaluator
  ParameterList model_evaluator_list = richards_list.get<Teuchos::ParameterList>("Richards model evaluator");  
  Amanzi::RichardsModelEvaluator RME(&problem, model_evaluator_list, problem.Map(), FS);

  // create the time stepping object
  Teuchos::RCP<Teuchos::ParameterList> bdf2_list(new Teuchos::ParameterList);
  *bdf2_list = richards_list.get<Teuchos::ParameterList>("Time integrator");  
  BDF2::Dae TS(RME, problem.Map());
  TS.setParameterList(bdf2_list);

  // create the initial condition
  Epetra_Vector u(problem.Map());

  Epetra_Vector *ucells = problem.CreateCellView(u);
  for (int c=0; c<ucells->MyLength(); c++) {
    const Point& xc = mesh->cell_centroid(c);
    //(*ucells)[c] = 101325.0 - 9800 * (xc[2] + 62.0);
    (*ucells)[c] = xc[2] * (xc[2] + 10.0);
  }

  Epetra_Vector *ufaces = problem.CreateFaceView(u);
  for (int f=0; f<ufaces->MyLength(); f++) {
    const Point& xf = mesh->face_centroid(f);
    (*ufaces)[f] = 101325.0 - 9800 * (xf[2] + 62.0);
    (*ufaces)[f] = xf[2] * (xf[2] + 10.0);
  }
  
  // initialize the state
  S->update_pressure(*problem.CreateCellView(u));
  S->set_time(0.0);

  // set intial and final time
  double t0 = 0.0;
  double t1 = 100.0;

  // compute the initial udot
  Epetra_Vector udot(problem.Map());
  problem.compute_udot(t0, u, udot);
  udot.PutScalar(0.0);

  // intialize the state of the time stepper
  TS.set_initial_state(t0, u, udot);
 
  int errc;
  double hnext, h = 1.0e-5;  // initial time step
  RME.update_precon(t0, u, h, errc);

  // set up output
  std::string cgns_filename = "out.cgns";

  Amanzi::CGNS::create_mesh_file(*mesh, cgns_filename);
  Amanzi::CGNS::open_data_file(cgns_filename);
  Amanzi::CGNS::create_timestep(0.0, 0, Amanzi::AmanziMesh::CELL);

  Amanzi::CGNS::write_field_data(*(S->get_pressure()), "pressure");
  Amanzi::CGNS::write_field_data(*(S->get_vertical_permeability()), "vertical permeability");
  Amanzi::CGNS::write_field_data(*(S->get_horizontal_permeability()), "horizontal permeability");
  

  // iterate
  int i = 0;
  double tlast = t0;
  do {    
    TS.bdf2_step(h, 0.0, u, hnext);
    TS.commit_solution(h, u);

    S->advance_time(h);    
    S->update_pressure(*problem.CreateCellView(u));  // update only the cell-based pressure

    TS.write_bdf2_stepping_statistics();

    h = hnext;
    i++;

    tlast=TS.most_recent_time();

    if ( i%50 == 1 ) { 
      Amanzi::CGNS::open_data_file(cgns_filename);
      Amanzi::CGNS::create_timestep(S->get_time(), i, Amanzi::AmanziMesh::CELL);
      Amanzi::CGNS::write_field_data(*(S->get_pressure()), "pressure");
      Amanzi::CGNS::write_field_data(*(S->get_vertical_permeability()), "vertical permeability");
      Amanzi::CGNS::write_field_data(*(S->get_horizontal_permeability()), "horizontal permeability");
    }
  } while (t1 >= tlast);

  Amanzi::CGNS::open_data_file(cgns_filename);
  Amanzi::CGNS::create_timestep(S->get_time(), i, Amanzi::AmanziMesh::CELL);
  Amanzi::CGNS::write_field_data(*(S->get_pressure()), "pressure");
  Amanzi::CGNS::write_field_data(*(S->get_vertical_permeability()), "vertical permeability");
  Amanzi::CGNS::write_field_data(*(S->get_horizontal_permeability()), "horizontal permeability");
  for (int k=0; k<100; k++) std::cout << k << " " << u[k] << std::endl;
 
  delete comm;
}
