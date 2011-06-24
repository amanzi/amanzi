#include <iostream>
#include <cstdlib>
#include <cmath>
#include "UnitTest++.h"
#include <vector>

#include "MSTK_types.h"

#include "Mesh_MSTK.hh"
#include "State.hpp"
#include "Transport_PK.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"



TEST(ADVANCE_WITH_MSTK) {
  using namespace Teuchos;
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziTransport;

  std::cout << "================ TEST ADVANCE WITH MSTK ===================" << endl;
  int num_components = 3;
  RCP<Mesh> mesh = rcp(new Mesh_MSTK("../mesh/mesh_mstk/test/hex_4x4x4_ss.exo", MPI_COMM_WORLD));

  State mpc_state(num_components, mesh);
  RCP<Transport_State> TS = rcp(new Transport_State(mpc_state));
  double u[3] = {1, 0, 0};

  TS->analytic_darcy_flux(u);
  TS->analytic_porosity();
  TS->analytic_water_saturation();
  TS->analytic_water_density();

  // initialize a transport process kernel from a transport state
  ParameterList parameter_list;
  string xmlFileName = "test/test_mstk.xml";

  updateParametersFromXmlFile(xmlFileName, &parameter_list);
  Transport_PK TPK(parameter_list, TS);

  // advance the state
  double  dT = TPK.calculate_transport_dT();
  TPK.advance(dT);

  // printing cell concentration 
  int i, k;
  double T = 0.0;
  RCP<Transport_State> TS_next = TPK.get_transport_state_next();
  RCP<Epetra_MultiVector> tcc = TS->get_total_component_concentration();
  RCP<Epetra_MultiVector> tcc_next = TS_next->get_total_component_concentration();

  for (i=0; i<50; i++) {
    dT = TPK.calculate_transport_dT();
    TPK.advance(dT);
    T += dT;

    if (i < 10) {
      printf( "T=%6.2f  C_0(x):", T );
      for( int k=0; k<9; k++ ) printf("%7.4f", (*tcc_next)[0][k]); std::cout << endl;
    }
    *tcc = *tcc_next;
  }

  // check that the final state is constant
  for( int k=0; k<4; k++ ) 
    CHECK_CLOSE( (*tcc_next)[0][k], 1.0, 1e-6 );
}
 


