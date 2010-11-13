#include <iostream>
#include <cstdlib>
#include <cmath>
#include "UnitTest++.h"
#include <vector>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "Mesh_maps_simple.hh"
#include "State.hpp"
#include "Transport_PK.hpp"



/* test constructor of transport PK */
TEST(CONSTRUCTOR) {

  using namespace std;
  using namespace Teuchos;

  cout << "================ TEST XML FILE ===================" << endl;
#ifdef HAVE_MPI
  Epetra_MpiComm    *comm = new Epetra_MpiComm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm *comm = new Epetra_SerialComm();
#endif

  /* create a MPC state with three component */
  int num_components = 3;
  RCP<Mesh_maps_simple>  mesh = rcp( new Mesh_maps_simple(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1, 2, 1, comm) ); 

  State mpc_state( num_components, mesh );

  /* create a transport state from the MPC state */
  RCP<Transport_State>  TS = rcp( new Transport_State(mpc_state) );

  /* initialize a transport process kernel from a transport state */
  ParameterList TPK_list;
  string xmlFileName = "test/test_transport.xml";

  updateParametersFromXmlFile( xmlFileName, &TPK_list );
  Transport_PK  TPK( TPK_list, TS );

  /* read the CFL number from the parameter list with a default of 1.0 */
  double cfl = TPK.get_cfl();
  CHECK( 0 < cfl && cfl <= 1.0 );
  cout << "CFL = " << cfl << endl;

  delete comm;
}




/* check that calculated faces and volumes are positive */
TEST(FACES_VOLUMES) {

  using namespace std;
  using namespace Teuchos;

  cout << "================ TEST FACES AND VOLUMES ===================" << endl;
#ifdef HAVE_MPI
  Epetra_MpiComm    *comm = new Epetra_MpiComm( MPI_COMM_WORLD);
#else
  Epetra_SerialComm *comm = new Epetra_SerialComm();
#endif

  /* create a MPC state with three component */
  int num_components = 3;
  RCP<Mesh_maps_simple>  mesh = rcp( new Mesh_maps_simple(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1, 2, 1, comm) ); 

  State mpc_state( num_components, mesh );

  /* create a transport state from the MPC state */
  RCP<Transport_State>  TS = rcp( new Transport_State(mpc_state) );

  /* initialize a transport process kernel from a transport state */
  ParameterList  parameter_list;
  string xmlFileName = "test/test_transport.xml";

  updateParametersFromXmlFile( xmlFileName, &parameter_list );
  Transport_PK  TPK( parameter_list, TS );

  /* printing face areas */
  int  f;
  double area;
  const Epetra_Map & face_map = mesh->face_map( true );

  cout << "Face areas: ";
  for( f=face_map.MinLID(); f<=face_map.MaxLID(); f++ ) { 
     area = TPK.get_face_area( f );
     cout << area << " ";
     CHECK_CLOSE( area, 0.75, 0.25 );
  }
  cout << endl;

  /* printing cell volumes */
  int  c;
  double volume;
  const Epetra_Map & cell_map = mesh->cell_map( true );

  cout << "Cell volumes: ";
  for( c=cell_map.MinLID(); c<=cell_map.MaxLID(); c++ ) { 
     volume = TPK.get_cell_volume( c );
     cout << volume << " ";
     CHECK_EQUAL( volume, 0.5 );
  }
  cout << endl;

  delete comm;
}
 



/* test advance routine */
TEST(ADVANCE_WITH_SIMPLE) {

  using namespace std;
  using namespace Teuchos;

  cout << "================ TEST ADVANCE ===================" << endl;
#ifdef HAVE_MPI
  Epetra_MpiComm     *comm = new Epetra_MpiComm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm  *comm = new Epetra_SerialComm();
#endif

  /* create a MPC state with three component */
  int num_components = 3;
  RCP<Mesh_maps_simple>  mesh = rcp( new Mesh_maps_simple(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 20, 2, 2, comm) ); 

  State mpc_state( num_components, mesh );

  /* create a transport state from the MPC state */
  RCP<Transport_State>  TS = rcp( new Transport_State(mpc_state) );

  /* initialize a transport process kernel from a transport state */
  ParameterList parameter_list;
  string xmlFileName = "test/test_transport.xml";

  updateParametersFromXmlFile( xmlFileName, &parameter_list );
  Transport_PK  TPK( parameter_list, TS );

  /* create analytic Darcy flux */
  double u[3] = {1, 0, 0};

  TS->analytic_darcy_flux( u );
  TS->analytic_porosity();
  TS->analytic_water_saturation();
  TS->analytic_water_density();

  /* advance the state */
  int  i, k;
  double  T = 0.0;
  RCP<Transport_State> TS_next = TPK.get_transport_state_next();

  RCP<Epetra_MultiVector> tcc      = TS->get_total_component_concentration();
  RCP<Epetra_MultiVector> tcc_next = TS_next->get_total_component_concentration();

  for( i=0; i<100; i++ ) {
     double dT = TPK.calculate_transport_dT();
     TPK.advance( dT );
     T += dT;

     if ( i < 10 ) {
        printf( "T=%6.1f  C_0(x):", T );
        for( int k=0; k<15; k++ ) printf("%7.4f", (*tcc_next)[0][k]); cout << endl;
     }

     for( int k=0; k<19; k++ ) 
        CHECK( ((*tcc_next)[0][k] - (*tcc_next)[0][k+1]) > -1e-15 );

     *tcc = *tcc_next;
  }

  /* check that the final state is constant */
  for( int k=0; k<20; k++ ) 
     CHECK_CLOSE( (*tcc_next)[0][k], 1.0, 1e-6 );

  delete comm;
}
 



/* convergence analisis */
TEST(CONVERGENCE_ANALYSIS) {

  using namespace std;
  using namespace Teuchos;

  cout << "================ TEST CONVERGENCE ANALISYS ===================" << endl;
  Epetra_SerialComm  *comm = new Epetra_SerialComm();

  /* create a MPC state with three component */
  int num_components = 1;
  RCP<Mesh_maps_simple>  mesh = rcp( new Mesh_maps_simple(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 40, 2, 2, comm) ); 

  State mpc_state( num_components, mesh );

  /* create a transport state from the MPC state */
  RCP<Transport_State>  TS = rcp( new Transport_State(mpc_state) );

  /* initialize a transport process kernel from a transport state */
  ParameterList parameter_list;
  string xmlFileName = "test/test_transport.xml";

  updateParametersFromXmlFile( xmlFileName, &parameter_list );
  Transport_PK  TPK( parameter_list, TS );

  /* create analytic Darcy flux */
  double u[3] = {1, 0, 0};

  TS->analytic_darcy_flux( u );
  TS->analytic_total_component_concentration();
  TS->analytic_porosity( 1.0 );
  TS->analytic_water_saturation( 1.0 );
  TS->analytic_water_density( 1.0 );

  /* advance the state */
  int  i, k, iter = 0;
  double  T = 0.0, T1 = 0.5;
  RCP<Transport_State> TS_next = TPK.get_transport_state_next();

  RCP<Epetra_MultiVector> tcc      = TS->get_total_component_concentration();
  RCP<Epetra_MultiVector> tcc_next = TS_next->get_total_component_concentration();

  while( T < T1 ) {
     double dT = min( TPK.calculate_transport_dT(), T1 - T );
     TPK.advance( dT );
     T += dT;

     *tcc = *tcc_next;
     iter++;
  }

  /* calculate L1 error */
  double  L1, L2;
  TS->error_total_component_concentration( T, TPK.get_cell_volume(), &L1, &L2 );
  cout << "L1 error = " << L1 << "  L2 error = " << L2 << "  dT = " << T1 / iter << endl;

  delete comm;
}
