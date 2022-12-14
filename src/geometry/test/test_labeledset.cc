/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Rao Garimella
*/

//
// Unit test to check if a labeled set region can be constructed correctly

#include <UnitTest++.h>

#include <iostream>

#include "mpi.h"

#include "AmanziComm.hh"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Array.hpp"

#include "../Point.hh"
#include "../Region.hh"
#include "../RegionLabeledSet.hh"
#include "../RegionFactory.hh"

TEST(LABELEDSET_REGION)
{
  auto ecomm = Amanzi::getDefaultComm();

  // read the parameter list from input file

  std::string infilename = "test/labeledsetregion.xml";
  Teuchos::ParameterXMLFileReader xmlreader(infilename);

  Teuchos::ParameterList reg_spec(xmlreader.getParameters());

  for (auto i = reg_spec.begin(); i != reg_spec.end(); i++) {
    const std::string reg_name = reg_spec.name(i);
    const unsigned int reg_id = 9959; // something arbitrary

    Teuchos::ParameterList reg_params = reg_spec.sublist(reg_name);

    // Create a Labeled Set Region
    Teuchos::RCP<const Amanzi::AmanziGeometry::Region> reg =
      Amanzi::AmanziGeometry::createRegion(reg_spec.name(i), reg_id, reg_params, *ecomm);

    // See if we retrieved the name and id correctly
    CHECK_EQUAL(reg->get_name(), reg_name);
    CHECK_EQUAL(reg->get_id(), reg_id);

    // Get the entity type and mesh file name directly from the XML
    CHECK_EQUAL(reg_spec.isSublist(reg_spec.name(i)), true);

    Teuchos::ParameterList::ConstIterator j = reg_params.begin();
    Teuchos::ParameterList labset_params = reg_params.sublist(reg_params.name(j));
    std::string in_entity_str = labset_params.get<std::string>("entity");


    // Make sure that the region type is a Labeled Set
    CHECK_EQUAL(reg->get_type(), Amanzi::AmanziGeometry::RegionType::LABELEDSET);

    // See if the min-max of the region were correctly retrieved
    Amanzi::AmanziGeometry::Point p, n;
    Teuchos::RCP<const Amanzi::AmanziGeometry::RegionLabeledSet> lsreg =
      Teuchos::rcp_dynamic_cast<const Amanzi::AmanziGeometry::RegionLabeledSet>(reg);

    // Did we get the entity string right?
    CHECK_EQUAL(in_entity_str, lsreg->entity_str());
  }


  // cannot test functionality without a mesh
}
