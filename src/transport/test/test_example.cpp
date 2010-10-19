#include <iostream>
#include "stdlib.h"
#include "math.h"
#include "UnitTest++.h"
#include <vector>

#include "Teuchos_RCP.hpp"

#include "../stk_mesh/Mesh_maps_stk.hh"
#include "../mpc/State.hpp"
#include "../Transport_PK.hpp"



TEST(TRANSPORT_GENERIC) {

  using namespace std;
  using namespace STK_mesh;

  /* create a state with 1 component */
  Teuchos::RCP<STK_mesh::Mesh_maps_stk>  mesh_amanzi;

  int number_components = 1;
  State global_state ( number_components, mesh_amanzi ) ;

  std::cout << "Hello" << std::endl;
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

