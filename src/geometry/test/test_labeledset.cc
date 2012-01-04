
//
// Unit test to check if a labeled set region can be constructed correctly
// Author: Rao Garimella
//

#include <UnitTest++.h>

#include <iostream>


#include "../Region.hh"
#include "../LabeledSetRegion.hh"
#include "../RegionFactory.hh"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Array.hpp"

#include "mpi.h"


TEST(LABELEDSET_REGION)
{


  // read the parameter list from input file

  std::string infilename = "test/labeledsetregion.xml";
  Teuchos::ParameterXMLFileReader xmlreader(infilename);

  Teuchos::ParameterList reg_spec(xmlreader.getParameters());

  for (Teuchos::ParameterList::ConstIterator i = reg_spec.begin(); 
       i != reg_spec.end(); i++) {

    const std::string reg_name = reg_spec.name(i);     
    const unsigned int reg_id = 9959;                   // something arbitrary

    Teuchos::ParameterList reg_params = reg_spec.sublist(reg_name);

    // Create a Labeled Set Region
  
    Amanzi::AmanziGeometry::RegionPtr reg = 
      Amanzi::AmanziGeometry::RegionFactory(reg_spec.name(i), reg_id, reg_params);
  
    // See if we retrieved the name and id correctly
  
    CHECK_EQUAL(reg->name(),reg_name);
    CHECK_EQUAL(reg->id(),reg_id);
  

    // Get the entity type and mesh file name directly from the XML
  

    CHECK_EQUAL(reg_spec.isSublist(reg_spec.name(i)),true);

  
    Teuchos::ParameterList::ConstIterator j = reg_params.begin();
    Teuchos::ParameterList labset_params = reg_params.sublist(reg_params.name(j));
    std::string in_entity_str = labset_params.get< std::string >("Entity");

    
    // Make sure that the region type is a Labeled Set

    CHECK_EQUAL(reg->type(),Amanzi::AmanziGeometry::LABELEDSET);
  
    // See if the min-max of the region were correctly retrieved
  
    Amanzi::AmanziGeometry::Point p, n;
      Amanzi::AmanziGeometry::LabeledSetRegionPtr lsreg =
      dynamic_cast<Amanzi::AmanziGeometry::LabeledSetRegionPtr> (reg);

    // Did we get the entity string right?

    CHECK_EQUAL(in_entity_str,lsreg->entity_str());

  }
}  



