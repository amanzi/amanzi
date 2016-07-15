#include "UnitTest++.h"

#include "Units.hh"

using namespace Amanzi;
using namespace Amanzi::Utils;

using namespace boost::units;

TEST(UNITS_LENGTH) 
{
  Units units("molar");
  double val(1.0);
  quantity<si::length> len = val * units.length_["cm"];
  std::cout << len << std::endl;
}


