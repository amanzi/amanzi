#include <iostream>
#include "stdlib.h"
#include "math.h"
#include "UnitTest++.h"
#include <vector>

TEST(EXAMPLE1) {

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



TEST(EXAMPLE2) {

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
TEST(EXAMPLE3) {

  using namespace std;

  vector<int> x(10);
  vector<int> y(10);

  x[1] = 1.0;
  y[1] = 2.0;

  CHECK_EQUAL(x[1],y[1]);
  CHECK_CLOSE(x[1],y[1],0.0001);
}
