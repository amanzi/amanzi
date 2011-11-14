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
#include "Point.hh"

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
TEST(DISPERSION) {
  using namespace Teuchos;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziTransport;
  using namespace Amanzi::AmanziGeometry;

  std::cout << "=== DISPERSION ===" << endl;
  Epetra_SerialComm  *comm = new Epetra_SerialComm();

  int nx = 20;
  RCP<Mesh> mesh = rcp(new Mesh_simple(0.0, 0.0, 0.0, 5.0, 1.0, 1.0, nx, 1, 1, comm)); 

  // create a MPC and Transport states with one component
  int num_components = 1;
  State mpc_state(num_components, mesh);
  RCP<Transport_State> TS = rcp(new Transport_State(mpc_state));

  Point u(1.0, 0.0, 0.0);
  TS->analytic_darcy_flux(u);
  TS->analytic_total_component_concentration(f_step);
  TS->analytic_porosity(1.0);
  TS->analytic_water_saturation(1.0);
  TS->analytic_water_density(1.0);

  // initialize a transport process kernel from a transport state
  ParameterList parameter_list;
  string xmlFileName = "test/test_dispersion.xml";

  updateParametersFromXmlFile(xmlFileName, &parameter_list);
  Transport_PK TPK(parameter_list, TS);
  TPK.print_statistics();
  TPK.verbosity_level = 0;

  // advance the state
  int i, k, iter = 0;
  double T = 0.0, T1 = 1.0;

  RCP<Transport_State> TS_next = TPK.get_transport_state_next();
  RCP<Epetra_MultiVector> tcc = TS->get_total_component_concentration();
  RCP<Epetra_MultiVector> tcc_next = TS_next->get_total_component_concentration();

  double dT, dT0;
  dT0 = TPK.calculate_transport_dT();

  while (T < T1) {
    dT = std::min(TPK.calculate_transport_dT(), T1 - T);
    dT = std::min(dT, dT0);
    for (int k=0; k<nx; k++) printf("%10.8f\n", (*tcc_next)[0][k]); 
    printf("\n");

    TPK.advance(dT);
    T += dT;
    //TPK.check_tracer_bounds(*tcc_next, 0, 0.0, 1.0, 1e-12);

    *tcc = *tcc_next;
    iter++;
  }
 
  delete comm;
}



