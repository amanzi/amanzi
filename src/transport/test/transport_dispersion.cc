#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>

#include "UnitTest++.h"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"

#include "MeshFactory.hh"
#include "MeshAudit.hh"
#include "Point.hh"

#include "State.hh"
#include "Transport_PK.hh"


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
TEST(DISPERSION) {
  using namespace Teuchos;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziTransport;
  using namespace Amanzi::AmanziGeometry;

  std::cout << "Test: dispersion" << endl;
#ifdef HAVE_MPI
  Epetra_MpiComm *comm = new Epetra_MpiComm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm *comm = new Epetra_SerialComm();
#endif

  // read parameter list
  ParameterList parameter_list;
  string xmlFileName = "test/transport_dispersion.xml";
  ParameterXMLFileReader xmlreader(xmlFileName);
  parameter_list = xmlreader.getParameters();

  // create an MSTK mesh framework
  ParameterList region_list = parameter_list.get<Teuchos::ParameterList>("Regions");
  GeometricModelPtr gm = new GeometricModel(3, region_list, comm);

  FrameworkPreference pref;
  pref.clear();
  pref.push_back(MSTK);

  MeshFactory factory(comm);
  factory.preference(pref);
  int nx = 20;
  RCP<Mesh> mesh = factory(0.0,0.0,0.0, 5.0,1.0,1.0, nx, 1, 1, gm); 

  // create a transport states with one component
  int num_components = 1;
  State mpc_state(num_components, 0, mesh);
  RCP<Transport_State> TS = rcp(new Transport_State(mpc_state));

  Point u(1.0, 0.0, 0.0);
  TS->AnalyticDarcyFlux(u);
  TS->AnalyticTotalComponentConcentration(f_step);
  TS->AnalyticPorosity(1.0);
  TS->AnalyticWaterSaturation(1.0);
  TS->AnalyticWaterDensity(1.0);

  // create transport PK  
  ParameterList transport_list = parameter_list.get<Teuchos::ParameterList>("Transport");
  Transport_PK TPK(transport_list, TS);
  TPK.InitPK();
  TPK.PrintStatistics();
  TPK.verbosity = TRANSPORT_VERBOSITY_NONE;

  // advance the state
  int i, k, iter = 0;
  double T = 0.0, T1 = 1.0;

  RCP<Transport_State> TS_next = TPK.transport_state_next();
  RCP<Epetra_MultiVector> tcc = TS->total_component_concentration();
  RCP<Epetra_MultiVector> tcc_next = TS_next->total_component_concentration();

  double dT, dT0;
  dT0 = TPK.CalculateTransportDt();

  while (T < T1) {
    dT = std::min(TPK.CalculateTransportDt(), T1 - T);
    dT = std::min(dT, dT0);
    // for (int k=0; k<nx; k++) printf("%10.8f\n", (*tcc_next)[0][k]); 
    // printf("\n");

    TPK.Advance(dT);
    T += dT;
    TPK.CheckTracerBounds(*tcc_next, 0, 0.0, 1.0, 1e-12);

    *tcc = *tcc_next;
    iter++;
  }
 
  delete comm;
}



