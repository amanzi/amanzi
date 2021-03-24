
//
// Unit test to check if a enumerated set region can be constructed correctly
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
#include "../RegionEnumerated.hh"
#include "../RegionFactory.hh"


TEST(ENUMERATEDSET_REGION)
{

  auto ecomm = Amanzi::getDefaultComm();

  // read the parameter list from input file

  std::string infilename = "test/enumeratedsetregion.xml";
  Teuchos::ParameterXMLFileReader xmlreader(infilename);

  Teuchos::ParameterList reg_spec(xmlreader.getParameters());

  for (Teuchos::ParameterList::ConstIterator i = reg_spec.begin(); 
       i != reg_spec.end(); i++) {

    const std::string reg_name = reg_spec.name(i);     
    const unsigned int reg_id = 9959;                   // something arbitrary

    Teuchos::ParameterList reg_params = reg_spec.sublist(reg_name);

    // Create a Enumerated Set Region
  
    Teuchos::RCP<const Amanzi::AmanziGeometry::Region> reg = 
      Amanzi::AmanziGeometry::createRegion(reg_spec.name(i), reg_id,
					   reg_params, *ecomm);
  
    // See if we retrieved the name and id correctly
    CHECK_EQUAL(reg->get_name(),reg_name);
    CHECK_EQUAL(reg->get_id(),reg_id);

    // Get the entity type and mesh file name directly from the XML
    CHECK_EQUAL(reg_spec.isSublist(reg_spec.name(i)),true);
  
    Teuchos::ParameterList::ConstIterator j = reg_params.begin();
    Teuchos::ParameterList labset_params = reg_params.sublist(reg_params.name(j));
    std::string in_entity_str = labset_params.get< std::string >("entity");
    
    // Make sure that the region type is a Enumerated Set
    CHECK_EQUAL(reg->get_type(),Amanzi::AmanziGeometry::ENUMERATED);
  
    // See if the min-max of the region were correctly retrieved
    Amanzi::AmanziGeometry::Point p, n;
    Teuchos::RCP<const Amanzi::AmanziGeometry::RegionEnumerated> lsreg =
      Teuchos::rcp_dynamic_cast<const Amanzi::AmanziGeometry::RegionEnumerated>(reg);

    // Did we get the entity string right?
    CHECK_EQUAL(in_entity_str,lsreg->entity_str());
  }

  // cannot test the functionality of the region without a mesh
}  



