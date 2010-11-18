#include <iostream>
#include <cstdlib>
#include <cmath>
#include "UnitTest++.h"
#include <vector>

#include "Mesh_maps_moab.hh"
#include "State.hpp"
#include "Transport_PK.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"


double f_step( double* x, double t ) { 
  if ( x[0] <= t ) return 1;
  return 0;
}


TEST(ADVANCE_WITH_MOAB_PARALLEL) {

  using namespace std;
  using namespace Teuchos;

  cout << "================ TEST PARALLEL MOAB MESH ===================" << endl;
  /* create a MPC state with one component */
  int num_components = 2;
  RCP<Mesh_maps_moab> mesh = rcp( new Mesh_maps_moab( "../moab_mesh/test/hex_4x4x4_ss_4P.h5m", MPI_COMM_WORLD ) );

  State mpc_state( num_components, mesh );

  /* create a transport state from the MPC state and populate it */
  RCP<Transport_State>  TS = rcp( new Transport_State( mpc_state ) );
  double u[3] = {1, 0, 0};

  TS->analytic_total_component_concentration( f_step );
  TS->analytic_porosity();
  TS->analytic_darcy_flux( u );
  TS->analytic_water_saturation();


  /* initialize a transport process kernel from a transport state */
  ParameterList parameter_list;
  string xmlFileName = "test/test_moab_parallel.xml";

  updateParametersFromXmlFile( xmlFileName, &parameter_list );
  Transport_PK  TPK( parameter_list, TS );

  /* advance the state */
  double  dT = TPK.calculate_transport_dT();
  TPK.advance( dT );

  /* printing cell concentration */
  int  iter, k;
  double  T = 0.0;
  RCP<Transport_State> TS_next = TPK.get_transport_state_next();

  RCP<Epetra_MultiVector> tcc      = TS->get_total_component_concentration();
  RCP<Epetra_MultiVector> tcc_next = TS_next->get_total_component_concentration();

  iter = 0;
  while( T < 1.0 ) {
     dT = TPK.calculate_transport_dT();
     TPK.advance( dT );
     T += dT;
     iter++;

     if ( iter < 10 && TPK.MyPID == 3 ) {
        printf( "T=%7.2f  C_0(x):", T );
        for( int k=0; k<10; k++ ) printf("%7.4f", (*tcc_next)[0][k]); cout << endl;
     }

     *tcc = *tcc_next;
  }

  /* check that the final state is constant */
  for( int k=0; k<4; k++ ) 
     CHECK_CLOSE( (*tcc_next)[0][k], 1.0, 1e-6 );
}
 
 


