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
  double z = (x[0]-1-t) / 2;
  return 2*z*z*z - 3*z*z + 1;
}


/* **************************************************************** */
TEST(CONVERGENCE_ANALYSIS_1ST) {
  using namespace Teuchos;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziTransport;
  using namespace Amanzi::AmanziGeometry;

  std::cout << "=== TEST CONVERGENCE ANALISYS 2ND ===" << endl;
  Epetra_SerialComm  *comm = new Epetra_SerialComm();

  for (int nx=10; nx<41; nx*=2 ) {
    RCP<Mesh> mesh = rcp(new Mesh_simple(0.0, 0.0, 0.0, 5.0, 1.0, 1.0, nx, 2, 2, comm)); 

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
    string xmlFileName = "test/test_high_order.xml";

    updateParametersFromXmlFile(xmlFileName, &parameter_list);
    Transport_PK TPK(parameter_list, TS);

    if (nx == 10) TPK.print_statistics();
    TPK.verbosity_level = 0;

    // advance the state
    int i, k, iter = 0;
    double T = 0.0, T1 = 2.0;

    RCP<Transport_State>    TS_next  = TPK.get_transport_state_next();
    RCP<Epetra_MultiVector> tcc      = TS->get_total_component_concentration();
    RCP<Epetra_MultiVector> tcc_next = TS_next->get_total_component_concentration();

    double dT, dT0;
    if (nx==10) dT0 = TPK.calculate_transport_dT();
    else dT0 /= 2;

    while (T < T1) {
      dT = std::min(TPK.calculate_transport_dT(), T1 - T);
      dT = std::min(dT, dT0);

      TPK.advance(dT);
      T += dT;
      if (TPK.internal_tests) {
        TPK.check_tracer_bounds(*tcc_next, 0, 0.0, 1.0, 1e-12);
      }

      *tcc = *tcc_next;
      iter++;
    }
    for (int k=0; k<nx; k++) cout << (*tcc_next)[0][k] << endl;

    double L1, L2;  // L1 and L2 errors
    TS->error_total_component_concentration(f_cubic, T, &L1, &L2);
    printf("nx=%3d  L1 error=%10.8f  L2 error=%10.8f  dT=%7.4f\n", nx, L1, L2, T1 / iter);
  }

  delete comm;
}



