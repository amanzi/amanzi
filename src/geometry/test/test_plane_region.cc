
//
// Unit test to check if a plane region can be constructed correctly
// Author: Rao Garimella
//

#include <UnitTest++.h>

#include <iostream>


#include "../Region.hh"
#include "../PlaneRegion.hh"
#include "../RegionFactory.hh"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Array.hpp"

#include "mpi.h"


TEST(PLANE_REGION)
{


  // read the parameter list from input file

  std::string infilename = "test/planeregion.xml";
  Teuchos::ParameterXMLFileReader xmlreader(infilename);

  Teuchos::ParameterList reg_spec(xmlreader.getParameters());


  Teuchos::ParameterList::ConstIterator i = reg_spec.begin();
  const std::string reg_name = reg_spec.name(i);     
  const unsigned int reg_id = 9959;                   // something arbitrary

  Teuchos::ParameterList reg_params = reg_spec.sublist(reg_name);

  // Create a rectangular region
  
  Amanzi::AmanziGeometry::RegionPtr reg = 
    Amanzi::AmanziGeometry::RegionFactory(reg_spec.name(i), reg_id, reg_params);
  
  // See if we retrieved the name and id correctly
  
  CHECK_EQUAL(reg->name(),reg_name);
  CHECK_EQUAL(reg->id(),reg_id);
  
  
  // Get the min-max bounds of the region from the XML specification
  
  Teuchos::Array<double> in_xyz, in_nrm;

  CHECK_EQUAL(reg_spec.isSublist(reg_spec.name(i)),true);

  Teuchos::ParameterList::ConstIterator j = reg_params.begin();
  Teuchos::ParameterList plane_params = reg_params.sublist(reg_params.name(j));
  in_xyz = plane_params.get< Teuchos::Array<double> >("Location");
  in_nrm = plane_params.get< Teuchos::Array<double> >("Direction");

  
 
  // Make sure that the region type is a Plane

  CHECK_EQUAL(reg->type(),Amanzi::AmanziGeometry::PLANE);
  
  // See if the min-max of the region were correctly retrieved
  
  Amanzi::AmanziGeometry::Point p, n;
  Amanzi::AmanziGeometry::PlaneRegionPtr plane =
    dynamic_cast<Amanzi::AmanziGeometry::PlaneRegionPtr> (reg);

  p = plane->point();
  n = plane->normal();

 
  CHECK_EQUAL(p.x(),in_xyz[0]);
  CHECK_EQUAL(p.y(),in_xyz[1]);
  CHECK_EQUAL(p.z(),in_xyz[2]);
  CHECK_EQUAL(n.x(),in_nrm[0]);
  CHECK_EQUAL(n.y(),in_nrm[1]);
  CHECK_EQUAL(n.z(),in_nrm[2]);
 
}  



