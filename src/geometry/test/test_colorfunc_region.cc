
//
// Unit test to check if a labeled set region can be constructed correctly
// Author: Rao Garimella
//


#include <iostream>
#include "mpi.h"

#include "UnitTest++.h"
#include "Teuchos_DefaultMpiComm.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Array.hpp"

#include "AmanziTypes.hh"

#include "../Point.hh"
#include "../Region.hh"
#include "../RegionColorFunction.hh"
#include "../RegionFactory.hh"

using namespace Amanzi;

TEST(COLORFUNCTION_REGION)
{
  auto ecomm = Comm_ptr_type(new Teuchos::MpiComm<int>(MPI_COMM_WORLD));

  // read the parameter list from input file
  std::string infilename = "test/colorfunc_region.xml";
  Teuchos::ParameterXMLFileReader xmlreader(infilename);

  Teuchos::ParameterList reg_spec(xmlreader.getParameters());

  for (auto i = reg_spec.begin(); i != reg_spec.end(); i++) {
    const std::string reg_name = reg_spec.name(i);     
    const unsigned int reg_id = 9959;                   // something arbitrary

    Teuchos::ParameterList reg_params = reg_spec.sublist(reg_name);

    // Create a Color Function Region
    Teuchos::RCP<const AmanziGeometry::Region> reg = 
      AmanziGeometry::createRegion(reg_spec.name(i), reg_id,
					   reg_params, ecomm);
  
    // See if we retrieved the name and id correctly
    CHECK_EQUAL(reg->name(),reg_name);
    CHECK_EQUAL(reg->id(),reg_id);
  
    // Make sure that the region type is an Indicator Function
    CHECK_EQUAL(reg->type(),AmanziGeometry::COLORFUNCTION);

    // Check if two known points are in the appropriate regions
    AmanziGeometry::Point p(3);

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



