/*
The transport component of the Amanzi code, serial unit tests.
License: BSD
Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>

#include "UnitTest++.h"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "Mesh_simple.hh"
#include "MeshAudit.hh"

#include "State.hpp"
#include "Transport_PK.hpp"


double f_step(double* x, double t) { 
  if (x[0] <= 1 + t) return 1;
  return 0;
}

double f_smooth(double* x, double t) { 
  return 0.5 - atan(50*(x[0]-5-t)) / M_PI;
}

double f_cubic(double* x, double t) { 
  if( x[0] < 1 + t ) return 1;
  if( x[0] > 3 + t ) return 0;
  double z = (*x-1-t) / 2;
  return 2*z*z*z - 3*z*z + 1;
}


/* **************************************************************** */
TEST(CONSTRUCTOR) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziTransport;

  std::cout << "=== TEST XML FILE ===" << endl;
#ifdef HAVE_MPI
  Epetra_MpiComm *comm = new Epetra_MpiComm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm *comm = new Epetra_SerialComm();
#endif

  int num_components = 3;
  RCP<Mesh> mesh = rcp(new Mesh_simple(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1, 2, 1, comm)); 
 
  MeshAudit audit(mesh);
  audit.Verify();

  State mpc_state(num_components, mesh);
  RCP<Transport_State> TS = rcp(new Transport_State(mpc_state));

  ParameterList TPK_list;
  string xmlFileName = "test/test_transport.xml";

  updateParametersFromXmlFile(xmlFileName, &TPK_list);
  Transport_PK TPK(TPK_list, TS);

  TPK.print_statistics();

  double cfl = TPK.get_cfl();
  CHECK(0 < cfl && cfl <= 1.0);
  std::cout << "CFL = " << cfl << endl;

  delete comm;
}


/* **************************************************************** */
TEST(FACES_VOLUMES) {
  using namespace Teuchos;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziTransport;

  std::cout << "=== TEST FACES AND VOLUMES ===" << endl;
#ifdef HAVE_MPI
  Epetra_MpiComm *comm = new Epetra_MpiComm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm *comm = new Epetra_SerialComm();
#endif

  int num_components = 3;
  RCP<Mesh> mesh = rcp(new Mesh_simple(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1, 2, 1, comm)); 

  State mpc_state(num_components, mesh);
  RCP<Transport_State> TS = rcp(new Transport_State(mpc_state));

  /* initialize a transport process kernel from a transport state */
  ParameterList  parameter_list;
  string xmlFileName = "test/test_transport.xml";

  updateParametersFromXmlFile(xmlFileName, &parameter_list);
  Transport_PK TPK(parameter_list, TS);

  TPK.print_statistics();

  /* printing face areas */
  int f;
  double area;
  const Epetra_Map& face_map = mesh->face_map(true);

  cout << "Face areas: ";
  for (f=face_map.MinLID(); f<=face_map.MaxLID(); f++) { 
    area = mesh->face_area(f);
    std::cout << area << " ";
    CHECK_CLOSE(area, 0.75, 0.25);
  }
  std::cout << endl;

  /* printing cell volumes */
  int  c;
  double volume;
  const Epetra_Map& cell_map = mesh->cell_map(true);

  std::cout << "Cell volumes: ";
  for (c=cell_map.MinLID(); c<=cell_map.MaxLID(); c++) { 
    volume = mesh->cell_volume(c);
    std::cout << volume << " ";
    CHECK_EQUAL(volume, 0.5);
  }
  std::cout << endl;

  delete comm;
}
 

/* **************************************************************** */
TEST(ADVANCE_WITH_SIMPLE) {
  using namespace Teuchos;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziTransport;
  using namespace Amanzi::AmanziGeometry;

  std::cout << "=== TEST ADVANCE ===" << endl;
#ifdef HAVE_MPI
  Epetra_MpiComm  *comm = new Epetra_MpiComm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm  *comm = new Epetra_SerialComm();
#endif

  /* create a MPC state with three component */
  int num_components = 3;
  RCP<Mesh> mesh = rcp(new Mesh_simple(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 20, 2, 2, comm)); 

  State mpc_state(num_components, mesh);

  /* create a transport state from the MPC state and populate it */
  RCP<Transport_State> TS = rcp(new Transport_State(mpc_state));
  Point u(1.0, 0.0, 0.0);

  TS->analytic_darcy_flux(u);
  TS->analytic_porosity();
  TS->analytic_water_saturation();
  TS->analytic_water_density();

  /* initialize a transport process kernel from a transport state */
  ParameterList parameter_list;
  string xmlFileName = "test/test_transport.xml";

  updateParametersFromXmlFile(xmlFileName, &parameter_list);
  Transport_PK TPK(parameter_list, TS);

  TPK.print_statistics();

  /* advance the state */
  int iter, k;
  double T = 0.0;
  RCP<Transport_State> TS_next = TPK.get_transport_state_next();

  RCP<Epetra_MultiVector> tcc      = TS->get_total_component_concentration();
  RCP<Epetra_MultiVector> tcc_next = TS_next->get_total_component_concentration();

  iter = 0;
  while (T < 1.0) {
    double dT = TPK.calculate_transport_dT();
    TPK.advance(dT);
    T += dT;
    iter++;

    if (iter < 10) {
      printf( "T=%6.2f  C_0(x):", T );
      for( int k=0; k<15; k++ ) printf("%7.4f", (*tcc_next)[0][k]); cout << endl;
    }

    for( int k=0; k<19; k++ ) 
      CHECK( ((*tcc_next)[0][k] - (*tcc_next)[0][k+1]) > -1e-15 );

    *tcc = *tcc_next;
  }

  /* check that the final state is constant */
  for (int k=0; k<20; k++) 
    CHECK_CLOSE((*tcc_next)[0][k], 1.0, 1e-6);

  delete comm;
}


/* **************************************************************** */
TEST(CONVERGENCE_ANALYSIS_DONOR) {
  using namespace Teuchos;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziTransport;
  using namespace Amanzi::AmanziGeometry;

  std::cout << "=== TEST CONVERGENCE ANALISYS: DONOR ===" << endl;
  Epetra_SerialComm  *comm = new Epetra_SerialComm();

  for (int nx=20; nx<321; nx*=2 ) {
    RCP<Mesh> mesh = rcp(new Mesh_simple(0.0, 0.0, 0.0, 5.0, 1.0, 1.0, nx, 1, 1, comm)); 

    // create a MPC and Transport states with one component
    int num_components = 1;
    State mpc_state(num_components, mesh);
    RCP<Transport_State> TS = rcp(new Transport_State(mpc_state));

    Point u(1.0, 0.0, 0.0);
    TS->analytic_darcy_flux(u);
    TS->analytic_total_component_concentration(f_cubic);
    TS->analytic_porosity(1.0);
    TS->analytic_water_saturation(1.0);
    TS->analytic_water_density(1.0);

    // initialize a transport process kernel from a transport state
    ParameterList parameter_list;
    string xmlFileName = "test/test_transport.xml";

    updateParametersFromXmlFile(xmlFileName, &parameter_list);
    Transport_PK TPK(parameter_list, TS);

    if (nx == 20) TPK.print_statistics();
    TPK.verbosity_level = 0;

    // advance the state
    int i, k, iter = 0;
    double T = 0.0, T1 = 1.0;

    RCP<Transport_State>    TS_next  = TPK.get_transport_state_next();
    RCP<Epetra_MultiVector> tcc      = TS->get_total_component_concentration();
    RCP<Epetra_MultiVector> tcc_next = TS_next->get_total_component_concentration();

    while (T < T1) {
      double dT = std::min(TPK.calculate_transport_dT(), T1 - T);
      TPK.advance(dT);
      T += dT;

      *tcc = *tcc_next;
      iter++;
    }

    // calculate L1 and L2 errors
    double L1, L2;
    TS->error_total_component_concentration(f_cubic, T, &L1, &L2);
    printf("nx=%3d  L1 error=%7.5f  L2 error=%7.5f  dT=%7.4f\n", nx, L1, L2, T1 / iter);
  }

  delete comm;
}



