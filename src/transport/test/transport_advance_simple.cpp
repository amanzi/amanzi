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


double f_step(const Amanzi::AmanziGeometry::Point& x, double t) { 
  if (x[0] <= 1 + t) return 1;
  return 0;
}

double f_smooth(const Amanzi::AmanziGeometry::Point& x, double t) { 
  return 0.5 - atan(50*(x[0]-5-t)) / M_PI;
}

double f_cubic(const Amanzi::AmanziGeometry::Point& x, double t) { 
  if( x[0] < 1 + t ) return 1;
  if( x[0] > 3 + t ) return 0;
  double z = (x[0]-1-t) / 2;
  return 2*z*z*z - 3*z*z + 1;
}


/* **************************************************************** */
TEST(ADVANCE_WITH_SIMPLE) {
  using namespace Teuchos;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziTransport;
  using namespace Amanzi::AmanziGeometry;

  std::cout << "Test: advance using simple mesh" << endl;
#ifdef HAVE_MPI
  Epetra_MpiComm  *comm = new Epetra_MpiComm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm  *comm = new Epetra_SerialComm();
#endif

  // read parameter list
  ParameterList parameter_list;
  string xmlFileName = "test/transport_advance_simple.xml";
  updateParametersFromXmlFile(xmlFileName, &parameter_list);

  // create an SIMPLE mesh framework 
  ParameterList region_list = parameter_list.get<Teuchos::ParameterList>("Regions");
  GeometricModelPtr gm = new GeometricModel(3, region_list, (Epetra_MpiComm *)comm);
  RCP<Mesh> mesh = rcp(new Mesh_simple(0.0,0.0,0.0, 1.0,1.0,1.0, 20, 20, 2, comm, gm)); 
  
  // create a transport state with one component 
  int num_components = 1;
  State mpc_state(num_components, 0, mesh);
  RCP<Transport_State> TS = rcp(new Transport_State(mpc_state));
 
  Point u(1.0, 0.0, 0.0);
  TS->AnalyticDarcyFlux(u);
  TS->AnalyticPorosity();
  TS->AnalyticWaterSaturation();
  TS->AnalyticWaterDensity();

  ParameterList transport_list = parameter_list.get<Teuchos::ParameterList>("Transport");
  Transport_PK TPK(transport_list, TS);
  TPK.InitPK();
  TPK.PrintStatistics();

  // advance the state
  int iter, k;
  double T = 0.0;
  RCP<Transport_State> TS_next = TPK.transport_state_next();
  RCP<Epetra_MultiVector> tcc = TS->total_component_concentration();
  RCP<Epetra_MultiVector> tcc_next = TS_next->total_component_concentration();

  iter = 0;
  while (T < 1.0) {
    double dT = TPK.CalculateTransportDt();
    TPK.Advance(dT);
    T += dT;
    iter++;

    if (iter < 10) {
      printf("T=%6.2f  C_0(x):", T);
      for (int k = 0; k < 15; k++) printf("%7.4f", (*tcc_next)[0][k]); cout << endl;
    }

    for(int k = 0; k < 19; k++) 
      CHECK(((*tcc_next)[0][k] - (*tcc_next)[0][k+1]) > -1e-15);

    *tcc = *tcc_next;
  }

  for (int k = 0; k < 20; k++)  // check that the final state is constant
    CHECK_CLOSE((*tcc_next)[0][k], 1.0, 1e-6);

  delete comm;
}




