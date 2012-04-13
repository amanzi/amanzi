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

#include "MSTK_types.h"
#include "Mesh_MSTK.hh"

#include "State.hpp"
#include "Flow_State.hpp"
#include "RichardsProblem.hpp"
#include "RichardsModelEvaluator.hpp"
#include "BDF2_Dae.hpp"

#include "gmv_mesh.hh"


/* **************************************************************** */
TEST(FLOW_1D_RICHARDS) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;

cout << "Test: Tensor Richards, a cube model" << endl;
#ifdef HAVE_MPI
  Epetra_MpiComm  *comm = new Epetra_MpiComm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm  *comm = new Epetra_SerialComm();
#endif

  /* read parameter list */
  ParameterList parameter_list;
  string xmlFileName = "test/flow_richards_tensor.xml";
  updateParametersFromXmlFile(xmlFileName, &parameter_list);

  // create an SIMPLE mesh framework 
  ParameterList region_list = parameter_list.get<Teuchos::ParameterList>("Regions");
  GeometricModelPtr gm = new GeometricModel(3, region_list, (Epetra_MpiComm *)comm);

  RCP<AmanziMesh::Mesh> mesh = rcp(new Mesh_MSTK(0.0,0.0,0.0, 1.0,1.0,1.0, 2, 3, 5, comm, gm)); 

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
  richards_list.set("gravity", -gravity[2]);  // We need to pass the absolute value, |g|.

  RichardsProblem problem(mesh, richards_list);
  problem.set_flow_state(FS);
  problem.set_absolute_permeability(FS->get_vertical_permeability(), FS->get_horizontal_permeability());

  // create the Richards Model Evaluator
  ParameterList model_evaluator_list = richards_list.get<Teuchos::ParameterList>("Richards model evaluator");  
  Amanzi::RichardsModelEvaluator RME(&problem, model_evaluator_list, problem.Map(), FS);
  RME.update_norm(1.0e-3, 1.0e-3);

  // create the time stepping object
  Teuchos::RCP<Teuchos::ParameterList> bdf2_list(new Teuchos::ParameterList);
  *bdf2_list = richards_list.get<Teuchos::ParameterList>("Time integrator");  
  BDF2::Dae TS(RME, problem.Map());
  TS.setParameterList(bdf2_list);

  // calculate the constant Darcy mass velocity
  double rho = FS->get_fluid_density();
  double mu = FS->get_fluid_viscosity();
  const Epetra_Vector& kh = FS->get_horizontal_permeability();
  const Epetra_Vector& kv = FS->get_vertical_permeability();

  Point g(0.0, 0.0, gravity[2]); 
  Point K(kh[0], kh[0], kv[0]);  // model the permeability tensor
  Point u0(1.0, 2.0, 3.0);
  Point v0(3);

  for (int i=0; i<3; i++) v0[i] = -u0[i] / K[i];
  v0 *= mu / rho;
  v0 += g * rho;
  cout << "rho=" << rho << " mu=" << mu << endl;
  cout << "K=" << K << " g=" << g << endl;
  cout << "v0=" << v0 << endl;

  // create the initial condition
  Epetra_Vector u(problem.Map());
  Epetra_Vector *ucells = problem.CreateCellView(u);
  for (int c=0; c<ucells->MyLength(); c++) {
    const Point& xc = mesh->cell_centroid(c);
    (*ucells)[c] = v0 * xc;
  }

  Epetra_Vector *ufaces = problem.CreateFaceView(u);
  for (int f=0; f<ufaces->MyLength(); f++) {
    const Point& xf = mesh->face_centroid(f);
    (*ufaces)[f] = v0 * xf;
  }
  
  // initialize the state
  S->update_pressure(*ucells);
  S->set_time(0.0);

  // set intial and final time
  double t0 = 0.0;
  double t1 = 1.0e+11;

  // compute the initial udot
  Epetra_Vector udot(problem.Map());
  problem.compute_udot(t0, u, udot);
  udot.PutScalar(0.0);

  // intialize the state of the time stepper
  TS.set_initial_state(t0, u, udot);
 
  int errc;
  double hnext, h = 1.0e+7;  // initial time step
  RME.update_precon(t0, u, h, errc);

  // iterate
  int i = 0;
  double tlast = t0;
  do {    
    TS.bdf2_step(h, 0.0, u, hnext);
    TS.commit_solution(h, u);

    S->advance_time(h);    
    S->update_pressure(*ucells);  // update only the cell-based pressure

    TS.write_bdf2_stepping_statistics();

    h = hnext;
    i++;

    tlast=TS.most_recent_time();
  } while (t1 >= tlast && i < 60);

  // Derive the Richards fluxes on faces
  double l1_error;
  Epetra_Vector darcy_flux(*ufaces);
  problem.DeriveDarcyFlux(u, darcy_flux, l1_error);
  std::cout << "L1 norm of the Richards flux discrepancy = " << l1_error << std::endl;

  // check accuracy
  double err_p = 0.0, err_f = 0.0;
  for (int c=0; c<(*ucells).MyLength(); c++) {
    const Point& xc = mesh->cell_centroid(c);
    double p_exact = v0 * xc;
    //std::cout << c << " p_num=" << (*ucells)[c] << " p_ex=" << p_exact << std::endl;
    err_p += pow((*ucells)[c] - p_exact, 2.0);
  }
  for (int f=0; f<(*ufaces).MyLength(); f++) {
    const Point& xf = mesh->face_centroid(f);
    const Point normal = mesh->face_normal(f);
  
    double p_exact = v0 * xf;
    double f_exact = u0 * normal / rho;
    //cout << f << " " << xf 
    //          << "  p_num=" << (*ufaces)[f] << " p_ex=" << p_exact 
    //          << "  flux_num=" << darcy_flux[f] << " f_ex=" << f_exact << endl;
    err_f += pow(darcy_flux[f] - f_exact, 2.0);
  }
  cout << "ERRs=" << sqrt(err_p) << " " << sqrt(err_f) << endl;

  CHECK(err_p < 1e-8);
  CHECK(err_f < 1e-8); 

  delete comm;
}
