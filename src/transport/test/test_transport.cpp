#include <iostream>
#include <cstdlib>
#include <cmath>
#include "UnitTest++.h"
#include <vector>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "simple_mesh/Mesh_maps_simple.hh"
#include "mpc/State.hpp"
#include "transport/Transport_PK.hpp"



TEST(INPUT_XML) {

  using namespace std;
  using namespace Teuchos;

  string xmlInFileName = "test/test_transport.xml";

  ParameterList TPK_list;
  updateParametersFromXmlFile(xmlInFileName, &TPK_list);

  /* read the CFL number from the parameter list with a default of 1.0 */
  int number_components;
  double cfl;

  cfl = TPK_list.get<double>( "CFL", 1.0 );
  CHECK( 0 < cfl && cfl <= 1.0 );

  number_components = TPK_list.get<int>( "number of components" );
  CHECK( number_components > 0 );

  cout << "CFL = " << cfl << endl;
  cout << "Total number of components = " << number_components << endl;

  /* read number of boundary consitions */ 
  ParameterList TPK_BC_list;
  int i, nBCs;

  TPK_BC_list = TPK_list.get<ParameterList>("Transport BCs");
  nBCs = TPK_BC_list.get<int>("number of BCs");
  CHECK( nBCs > 0 );

  for( i=0; i<nBCs; i++ ) {
     char bc_char_name[10];
    
     sprintf(bc_char_name, "BC %d", i);

     string bc_name(bc_char_name);
     if ( ! TPK_BC_list.isSublist(bc_name) ) throw exception();

     ParameterList bc_list = TPK_BC_list.sublist(bc_name);

     int     ssid, ntcc;
     string  type;
     double  value;

     ssid  = bc_list.get<int>("Side set ID");
     ntcc  = bc_list.get<int>("number of components");
     type  = bc_list.get<string>("Type");
     value = bc_list.get<double>("Component 1");

     cout << "Boundary Condition: " << bc_name << endl;
     cout << "  side set id = " << ssid << endl;
     cout << "  type        = " << type << endl;
     cout << "  # componets = " << ntcc << endl;
     cout << "  component 2 = " << value << endl;
  }
  cout << "==================================================================" << endl << endl;
}



TEST(INIT_PROCESS_KERNEL) {

  using namespace std;
  using namespace Teuchos;

#ifdef HAVE_MPI
  Epetra_MpiComm    *comm = new Epetra_MpiComm( MPI_COMM_WORLD);
#else
  Epetra_SerialComm *comm = new Epetra_SerialComm();
#endif

  /* create a MPC state with one component */
  RCP<Mesh_maps_simple>  mesh = rcp( new Mesh_maps_simple(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1, 2, 1, comm) ); 

  int number_components = 3;
  State mpc_state( number_components, mesh );

  /*create a transport state from the MPC state */
  RCP<Transport_State>  TS = rcp( new Transport_State( mpc_state ) );

  /* read the trasport parameter list */
  string xmlInFileName = "test/test_transport.xml";

  ParameterList parameter_list;
  updateParametersFromXmlFile(xmlInFileName, &parameter_list);

  /* initialize a transport process kernel from a transport state */
  Transport_PK  TPK(parameter_list, TS);

  /* the actual test is to print the Darcy velocity */
  double u[3] = {1, 2, 3};
  TS->analytic_darcy_flux( u );

  RCP<Transport_State>  TS_pointer;
  RCP<const Epetra_Vector>  darcy_flux;

  TS_pointer = TPK.get_transport_state();
  darcy_flux = TS_pointer->get_darcy_flux();

  cout << *(darcy_flux) << endl;
  cout << "==================================================================" << endl << endl;

  delete comm;
}
 


TEST(FACES_VOLUMES) {

  using namespace std;
  using namespace Teuchos;

#ifdef HAVE_MPI
  Epetra_MpiComm    *comm = new Epetra_MpiComm( MPI_COMM_WORLD);
#else
  Epetra_SerialComm *comm = new Epetra_SerialComm();
#endif

  /* create a MPC state with one component */
  RCP<Mesh_maps_simple>  mesh = rcp( new Mesh_maps_simple(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1, 2, 1, comm) ); 

  int num_components = 3;
  State mpc_state( num_components, mesh );

  /* create a transport state from the MPC state */
  RCP<Transport_State>  TS = rcp( new Transport_State(mpc_state) );

  /* read the trasport parameter list */
  string xmlInFileName = "test/test_transport.xml";

  ParameterList parameter_list;
  updateParametersFromXmlFile(xmlInFileName, &parameter_list);

  /* initialize a transport process kernel from a transport state */
  Transport_PK  TPK(parameter_list, TS);

  /* printing face areas */
  int  f;
  double area;
  Epetra_Map face_map = mesh->face_map(false);

  for( f=face_map.MinLID(); f<=face_map.MaxLID(); f++ ) { 
     area = TPK.get_face_area( f );
     cout << "face id = " << f << "  area = " << area << endl;
     CHECK_CLOSE(area, 0.75, 0.25);
  }

  /* printing cell volumes */
  int  c;
  double volume;
  Epetra_Map cell_map = mesh->cell_map(false);

  for( c=cell_map.MinLID(); c<=cell_map.MaxLID(); c++ ) { 
     volume = TPK.get_cell_volume( c );
     cout << "cell id = " << c << "  volume = " << volume << endl;
     CHECK_EQUAL(area, 0.5);
  }
  cout << "==================================================================" << endl << endl;

  delete comm;
}
 


TEST(ADVANCE_WITH_SIMPLE) {

  using namespace std;
  using namespace Teuchos;

#ifdef HAVE_MPI
  Epetra_MpiComm    *comm = new Epetra_MpiComm( MPI_COMM_WORLD );
#else
  Epetra_SerialComm *comm = new Epetra_SerialComm();
#endif

  /* create a MPC state with one component */
  RCP<Mesh_maps_simple>  mesh = rcp( new Mesh_maps_simple(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 20, 2, 2, comm) ); 

  int num_components = 3;
  State mpc_state( num_components, mesh );

  /* create a transport state from the MPC state */
  RCP<Transport_State>  TS = rcp( new Transport_State(mpc_state) );

  /* read the trasport parameter list */
  string xmlInFileName = "test/test_transport.xml";

  ParameterList parameter_list;
  updateParametersFromXmlFile(xmlInFileName, &parameter_list);

  /* initialize a transport process kernel from a transport state */
  Transport_PK  TPK(parameter_list, TS);

  /* create analytic Darcy flux */
  double u[3] = {1, 0, 0};

  TS->analytic_total_component_concentration();
  TS->analytic_porosity();
  TS->analytic_darcy_flux( u );
  TS->analytic_water_saturation();
  TS->analytic_water_density();

  /* advance the state */
  TPK.advance();

  /* printing cell concentration */
  int  i, k;
  double  dT, T;
  RCP<Transport_State> TS_next = TPK.get_transport_state_next();

  RCP<Epetra_MultiVector> tcc      = TS->get_total_component_concentration();
  RCP<Epetra_MultiVector> tcc_next = TS_next->get_total_component_concentration();

  T = 0.0;
  cout << "Original state: T=" << T << endl;
  for( i=0; i<3; i++ ) {
     for( int k=0; k<20; k++ ) printf("%7.4f", (*tcc)[i][k]); 
     cout << endl;
  }
  dT = TPK.get_transport_dT();
  T += dT;

  cout << "New state: T=" << T << endl;
  for( i=0; i<3; i++ ) {
     for( int k=0; k<20; k++ ) printf("%7.4f", (*tcc_next)[i][k]); 
     cout << endl;
  }
 
  cout << "Dynamics of component 0 (3D simulaiton of 1D transport)" << endl;
  for( int i=0; i<20; i++ ) {
     *tcc = *tcc_next;

     TPK.advance();

     dT = TPK.get_transport_dT();
     T += dT;

     printf("T=%6.1f  C_0(x):", T);
     for( int k=0; k<20; k++ ) printf("%7.4f", (*tcc_next)[0][k]); cout << endl;
  }

  cout << "==================================================================" << endl << endl;

  delete comm;
}
 


