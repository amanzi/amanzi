#include "UnitTest++.h"

#include "ExprTK.hh"

using namespace Amanzi;
using namespace Amanzi::Utils;

TEST(EXPRTK) 
{
  std::string formula;
  std::vector<double> txyz;

  {
    ExprTK exprtk;
    txyz = { 0.5 };
    CHECK(exprtk.Initialize(1, "t"));
    CHECK_CLOSE(0.5, exprtk(txyz), 1e-8);
  }

  {
    ExprTK exprtk;
    txyz = { 0.5, 2.0 };
    CHECK(exprtk.Initialize(2, "t + x * x"));
    CHECK_CLOSE(4.5, exprtk(txyz), 1e-8);
  }

  {
    ExprTK exprtk;
    txyz = { 0.5, 2.0, 3.0 };
    CHECK(exprtk.Initialize(3, "t + x * x + sin(y)"));
    CHECK_CLOSE(4.641120008059867, exprtk(txyz), 1e-8);
  }

  {
    ExprTK exprtk;
    txyz = { 0.5, 2.0, 3.0, 4.0 };
    CHECK(exprtk.Initialize(4, "t + x * x + sin(y) + z / y^3"));
    CHECK_CLOSE(4.789268156208015, exprtk(txyz), 1e-8);
  }
}

