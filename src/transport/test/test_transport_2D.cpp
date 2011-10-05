/*
The transport component of the Amanzi code, serial unit tests.
License: BSD
Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>

#include "UnitTest++.h"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "MSTK_types.h"
#include "Mesh_MSTK.hh"
#include "MeshAudit.hh"

#include "State.hpp"
#include "Transport_PK.hpp"


/* **************************************************************** */
TEST(ADVANCE_WITH_2D_MESH) {
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
  RCP<Mesh> mesh = rcp(new Mesh_MSTK("test/rect2D_ss.exo", MPI_COMM_WORLD, 2));

  State mpc_state(num_components, mesh);

  /* create a transport state from the MPC state and populate it */
  RCP<Transport_State> TS = rcp(new Transport_State(mpc_state));
  Point u(1.0, 0.0);

  TS->analytic_darcy_flux(u);
  TS->analytic_porosity();
  TS->analytic_water_saturation();
  TS->analytic_water_density();

  /* initialize a transport process kernel from a transport state */
  ParameterList parameter_list;
  string xmlFileName = "test/test_transport_2D.xml";

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
      for( int k=0; k<9; k++ ) {
        int k1 = 9 - k;  // reflects cell numbering
        printf("%7.4f", (*tcc_next)[0][k1]); 
      }
      printf("\n");
    }

    for( int k=0; k<8; k++ ) {
      CHECK( ((*tcc_next)[0][k+1] - (*tcc_next)[0][k]) > -1e-15 );
    }

    *tcc = *tcc_next;
  }

  /* check that the final state is constant */
  for (int k=0; k<10; k++) 
    CHECK_CLOSE((*tcc_next)[0][k], 1.0, 1e-6);

  delete comm;
}





