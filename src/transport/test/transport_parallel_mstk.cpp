#include <iostream>
#include <cstdlib>
#include <cmath>
#include "UnitTest++.h"
#include <vector>

#include "Mesh_MSTK.hh"
#include "State.hpp"
#include "Transport_PK.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"


double f_step(const Amanzi::AmanziGeometry::Point& x, double t ) { 
  if ( x[0] <= t ) return 1;
  return 0;
}


TEST(ADVANCE_WITH_MSTK_PARALLEL) {
  using namespace std;
  using namespace Teuchos;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziTransport;
  using namespace Amanzi::AmanziGeometry;

  cout << "Test: advance with MSTK in parallel" << endl;
  // read parameter list
  ParameterList parameter_list;
  string xmlFileName = "test/transport_parallel_mstk.xml";
  updateParametersFromXmlFile(xmlFileName, &parameter_list);

  // create an MSTK mesh framework 
  ParameterList region_list = parameter_list.get<Teuchos::ParameterList>("Regions");
  GeometricModelPtr gm = new GeometricModel(3, region_list);
  RCP<Mesh> mesh = rcp(new Mesh_MSTK("../mesh/mesh_mstk/test/hex_4x4x4_ss.exo", MPI_COMM_WORLD, 3, gm));
   
  // create a transport state with two component 
  int num_components = 2;
  State mpc_state(num_components, mesh);
  RCP<Transport_State> TS = rcp(new Transport_State(mpc_state));

  Point u(1.0, 0.0, 0.0);
  TS->analytic_total_component_concentration(f_step);
  TS->analytic_porosity();
  TS->analytic_darcy_flux(u);
  TS->analytic_water_saturation();

  ParameterList transport_list =  parameter_list.get<Teuchos::ParameterList>("Transport");
  Transport_PK TPK(transport_list, TS);
  TPK.set_standalone_mode(true);
  TPK.print_statistics();

  // advance the state
  double dT = TPK.calculate_transport_dT();
  TPK.advance(dT);

  // printing cell concentration
  int  iter, k;
  double  T = 0.0;
  RCP<Transport_State> TS_next = TPK.get_transport_state_next();

  RCP<Epetra_MultiVector> tcc = TS->get_total_component_concentration();
  RCP<Epetra_MultiVector> tcc_next = TS_next->get_total_component_concentration();

  iter = 0;
  while(T < 1.0) {
    dT = TPK.calculate_transport_dT();
    TPK.advance(dT);
    T += dT;
    iter++;

    if (iter < 10 && TPK.MyPID == 3) {
      printf("T=%7.2f  C_0(x):", T);
      for (int k=0; k<2; k++) printf("%7.4f", (*tcc_next)[0][k]); cout << endl;
    }
    *tcc = *tcc_next;
  }

  //for (int k=0; k<12; k++) 
  //  CHECK_CLOSE((*tcc_next)[0][k], 1.0, 1e-6);
}
 
 


