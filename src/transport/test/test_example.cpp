#include <iostream>
#include "stdlib.h"
#include "math.h"
#include "UnitTest++.h"
#include <vector>

#include "Teuchos_RCP.hpp"

#include "../simple_mesh/Mesh_maps_simple.hh"
#include "../mpc/State.hpp"
#include "../Transport_PK.hpp"



TEST(TRANSPORT_PK_INIT) {

  using namespace std;
  using namespace Teuchos;

#ifdef HAVE_MPI
  Epetra_MpiComm    *comm = new Epetra_MpiComm( MPI_COMM_WORLD);
#else
  Epetra_SerialComm *comm = new Epetra_SerialComm();
#endif

  /* create MPC state with 1 component */
  RCP<Mesh_maps_simple>  mesh_amanzi = rcp( new Mesh_maps_simple(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1, 2, 1, comm) ); 

  int number_components = 1;
  State mpc_state( number_components, mesh_amanzi );

  /*create transport state from MPC state */
  RCP<Transport_State>  TS = rcp( new Transport_State(mpc_state) );

  /* initialize a transport process kernel from a transport state */
  Transport_PK  TPK(TS);


  /* the actual test is to print the Darcy velocity */
  RCP<Transport_State>  TS_pointer;
  RCP<Epetra_MultiVector>  darcy_flux_pointer;

  TS_pointer = TPK.get_transport_state();
  darcy_flux_pointer = TS_pointer->get_total_component_concentration();

  cout << *(darcy_flux_pointer) << endl;
}
 


TEST(EXAMPLE2) {

  using namespace std;

  vector<int> x(10);
  vector<int> y(10);

  for (int it = 0; it < 10; it++)
    {
      x[it] = 0.0;
    }

  for (int it = 0; it < 10; it++)
    {
      x[it] = 0.0;
    }
  

  CHECK_ARRAY_EQUAL(x,y,10);
  CHECK_ARRAY_CLOSE(x,y,10,0.00001);
}



TEST(EXAMPLE3) {

  using namespace std;

  vector<int> x(10);
  vector<int> y(10);

  x[1] = 1.0;
  y[1] = 1.0;
  y[2] = 2.0;

  CHECK_EQUAL(x[1],y[1]);
  CHECK_CLOSE(x[1],y[1],0.0001);
}



// this test will obviously fail
/*
TEST(EXAMPLE4) {

  using namespace std;

  vector<int> x(10);
  vector<int> y(10);

  x[1] = 1.0;
  y[1] = 2.0;

  CHECK_EQUAL(x[1],y[1]);
  CHECK_CLOSE(x[1],y[1],0.0001);
}
*/

