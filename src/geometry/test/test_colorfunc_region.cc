
//
// Unit test to check if a labeled set region can be constructed correctly
// Author: Rao Garimella
//

#include <UnitTest++.h>

#include <iostream>


#include "../Region.hh"
#include "../ColorFunctionRegion.hh"
#include "../RegionFactory.hh"

#include "Epetra_MpiComm.h"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Array.hpp"

#include "mpi.h"


TEST(COLORFUNCTION_REGION)
{

  Epetra_MpiComm ecomm(MPI_COMM_WORLD);


  // read the parameter list from input file

  std::string infilename = "test/colorfunc_region.xml";
  Teuchos::ParameterXMLFileReader xmlreader(infilename);

  Teuchos::ParameterList reg_spec(xmlreader.getParameters());

  for (Teuchos::ParameterList::ConstIterator i = reg_spec.begin(); 
       i != reg_spec.end(); i++) {

    const std::string reg_name = reg_spec.name(i);     
    const unsigned int reg_id = 9959;                   // something arbitrary

    Teuchos::ParameterList reg_params = reg_spec.sublist(reg_name);

    // Create a Labeled Set Region
  
    Amanzi::AmanziGeometry::RegionPtr reg = 
      Amanzi::AmanziGeometry::RegionFactory(reg_spec.name(i), reg_id, reg_params, &ecomm);
  
    // See if we retrieved the name and id correctly
  
    CHECK_EQUAL(reg->name(),reg_name);
    CHECK_EQUAL(reg->id(),reg_id);
  
    // Make sure that the region type is an Indicator Function

    CHECK_EQUAL(reg->type(),Amanzi::AmanziGeometry::COLORFUNCTION);

    // Check if two known points are in the appropriate regions

    Amanzi::AmanziGeometry::Point p(3);

    if (reg_name == "Top") {
      double xyz[3] = {0.5,0.5,0.75};
      p.set(xyz);

      CHECK_EQUAL(reg->inside(p),true);
    }
    else if (reg_name == "Bottom") {
      double xyz[3] = {0.5,0.5,0.25};
      p.set(xyz);

      CHECK_EQUAL(reg->inside(p),true);
    }
  }
}  



