#include <UnitTest++.h>

#include <iostream>
#include <vector>

#include "../Region.hh"
#include "../RegionCylinder.hh"
#include "../RegionFactory.hh"

#include "AmanziComm.hh"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Array.hpp"

TEST(CYLINDER_REGION)
{
  auto ecomm = Amanzi::getDefaultComm();

  // read the parameter list from input file
  std::string infilename = "test/cylinderregion.xml";
  Teuchos::ParameterXMLFileReader xmlreader(infilename);
  Teuchos::ParameterList reg_spec(xmlreader.getParameters());

  const std::string reg_name = reg_spec.name(reg_spec.begin());     
  const unsigned int reg_id = 9959;  // something arbitrary

  Teuchos::ParameterList reg_params = reg_spec.sublist(reg_name);
    
  // Create a rectangular region
  auto reg = Amanzi::AmanziGeometry::createRegion(reg_name, reg_id, reg_params, *ecomm);
  
  // See if we retrieved the name, id, and type correctly
  CHECK_EQUAL(reg->get_name(), reg_name);
  CHECK_EQUAL(reg->get_id(), reg_id);
  CHECK_EQUAL(reg_spec.isSublist(reg_spec.name(reg_spec.begin())), true);
  CHECK_EQUAL(reg->get_type(), Amanzi::AmanziGeometry::RegionType::CYLINDER);

  // Make sure that the region dimension is 3
  CHECK_EQUAL(reg->get_manifold_dimension(), 3);
  
  // test the functionality of the region
  std::vector<Amanzi::AmanziGeometry::Point> pin;
  pin.push_back(Amanzi::AmanziGeometry::Point(0.0, 0.0, 0.0));
  pin.push_back(Amanzi::AmanziGeometry::Point(0.0, 2.0, 1.0));
  pin.push_back(Amanzi::AmanziGeometry::Point(0.0,-2.0, 1.0));
  pin.push_back(Amanzi::AmanziGeometry::Point(0.5, 0.0, 2.0));

  for (auto p = pin.begin(); p != pin.end(); ++p) {
    CHECK(reg->inside(*p));
  }

  std::vector<Amanzi::AmanziGeometry::Point> pout;
  pin.push_back(Amanzi::AmanziGeometry::Point(9.9, 8.0, 0.0));
  pin.push_back(Amanzi::AmanziGeometry::Point(1.0, 2.0,-1.0));

  for (auto p = pout.begin(); p != pout.end(); ++p) {
    CHECK(!reg->inside(*p));
  }
}


