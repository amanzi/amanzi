
//
// Unit test to check if a labeled set region can be constructed correctly
// Author: Rao Garimella
//


#include <iostream>
#include "mpi.h"

#include "UnitTest++.h"
#include "AmanziComm.hh"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Array.hpp"

#include "../Point.hh"
#include "../Region.hh"
#include "../RegionFunctionColor.hh"
#include "../RegionFactory.hh"


TEST(COLORFUNCTION_REGION)
{
  auto ecomm = Amanzi::getDefaultComm();

  // read the parameter list from input file
  std::string infilename = "test/colorfunc_region.xml";
  Teuchos::ParameterXMLFileReader xmlreader(infilename);

  Teuchos::ParameterList reg_spec(xmlreader.getParameters());

  for (auto i = reg_spec.begin(); i != reg_spec.end(); i++) {
    const std::string reg_name = reg_spec.name(i);     
    const unsigned int reg_id = 9959;                   // something arbitrary

    Teuchos::ParameterList reg_params = reg_spec.sublist(reg_name);

    // Create a Color Function Region
    Teuchos::RCP<const Amanzi::AmanziGeometry::Region> reg = 
      Amanzi::AmanziGeometry::createRegion(reg_spec.name(i), reg_id,
					   reg_params, *ecomm);
  
    // See if we retrieved the name and id correctly
    CHECK_EQUAL(reg->get_name(),reg_name);
    CHECK_EQUAL(reg->get_id(),reg_id);
  
    // Make sure that the region type is an Indicator Function
    CHECK_EQUAL(reg->get_type(),Amanzi::AmanziGeometry::COLORFUNCTION);

    // Check if two known points are in the appropriate regions
    Amanzi::AmanziGeometry::Point p(3);

    if (reg_name == "Top") {
      double xyz[3] = {0.5,0.5,0.75};
      p.set(xyz);

      CHECK_EQUAL(reg->inside(p),true);

    } else if (reg_name == "Bottom") {
      double xyz[3] = {0.5,0.5,0.25};
      p.set(xyz);

      CHECK_EQUAL(reg->inside(p),true);
    }
  }
}  



