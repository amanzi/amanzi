#include <UnitTest++.h>

#include <iostream>


#include "../Region.hh"
#include "../RegionBox.hh"
#include "../RegionFactory.hh"

#include "Epetra_MpiComm.h"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Array.hpp"

#include "mpi.h"

TEST(BOX_REGION_2D)
{
  Epetra_MpiComm ecomm(MPI_COMM_WORLD);

  // read the parameter list from input file
  std::string infilename = "test/boxregion_2D.xml";
  Teuchos::ParameterXMLFileReader xmlreader(infilename);

  Teuchos::ParameterList reg_spec(xmlreader.getParameters());

  Teuchos::ParameterList::ConstIterator i = reg_spec.begin();
  const std::string reg_name = reg_spec.name(i);     
  const unsigned int reg_id = 9959;                   // something arbitrary

  Teuchos::ParameterList reg_params = reg_spec.sublist(reg_name);
    
  // Create a rectangular region
  Teuchos::RCP<Amanzi::AmanziGeometry::Region> reg = 
    Amanzi::AmanziGeometry::createRegion(reg_name, reg_id, reg_params, &ecomm);
  
  // See if we retrieved the name and id correctly
  CHECK_EQUAL(reg->name(),reg_name);
  CHECK_EQUAL(reg->id(),reg_id);
  CHECK_EQUAL(reg_spec.isSublist(reg_spec.name(i)),true);
  
  // Get the min-max bounds of the region from the XML specification
  Teuchos::Array<double> in_min_xyz, in_max_xyz;

  Teuchos::ParameterList::ConstIterator j = reg_params.begin();
  Teuchos::ParameterList box_params = reg_params.sublist(reg_params.name(j));
  in_min_xyz = box_params.get< Teuchos::Array<double> >("Low Coordinate");
  in_max_xyz = box_params.get< Teuchos::Array<double> >("High Coordinate");
 
  // Make sure that the region type is a BOX
  CHECK_EQUAL(reg->type(),Amanzi::AmanziGeometry::BOX);

  // Make sure that the region dimension is 2
  CHECK_EQUAL(reg->topological_dimension(),2);
  
  // See if the min-max of the region were correctly retrieved
  Amanzi::AmanziGeometry::Point pmin, pmax;
  Teuchos::RCP<const Amanzi::AmanziGeometry::RegionBox> rect =
    Teuchos::rcp_dynamic_cast<const Amanzi::AmanziGeometry::RegionBox>(reg);

  // Make sure we got back 2D points
  rect->corners(&pmin,&pmax);
  CHECK_EQUAL(pmin.dim(),2);
  CHECK_EQUAL(pmax.dim(),2);

  // Compare coordinates read from XML file and retrieved from region
  CHECK_EQUAL(pmin.x(),in_min_xyz[0]);
  CHECK_EQUAL(pmin.y(),in_min_xyz[1]);
  CHECK_EQUAL(pmax.x(),in_max_xyz[0]);
  CHECK_EQUAL(pmax.y(),in_max_xyz[1]);

  // test the functionality of the region
  std::vector<Amanzi::AmanziGeometry::Point> pin;
  pin.push_back(Amanzi::AmanziGeometry::Point(9.,8.));
  pin.push_back(Amanzi::AmanziGeometry::Point(11.,8.));
  pin.push_back(Amanzi::AmanziGeometry::Point(11.,1.5));
  pin.push_back(Amanzi::AmanziGeometry::Point(9.,1.5));
  pin.push_back(Amanzi::AmanziGeometry::Point(10., 5.));

  for (std::vector<Amanzi::AmanziGeometry::Point>::iterator p=pin.begin();
       p!=pin.end(); ++p) {
    CHECK(reg->inside(*p));
  }

  std::vector<Amanzi::AmanziGeometry::Point> pout;
  pin.push_back(Amanzi::AmanziGeometry::Point(9.9,8.));
  pin.push_back(Amanzi::AmanziGeometry::Point(11,7.9));

  for (std::vector<Amanzi::AmanziGeometry::Point>::iterator p=pout.begin();
       p!=pout.end(); ++p) {
    CHECK(!reg->inside(*p));
  }
}


TEST(BOX_REGION_3D)
{

  Epetra_MpiComm ecomm(MPI_COMM_WORLD);

  // read the parameter list from input file

  std::string infilename = "test/boxregion_3D.xml";
  Teuchos::ParameterXMLFileReader xmlreader(infilename);

  Teuchos::ParameterList reg_spec(xmlreader.getParameters());


  Teuchos::ParameterList::ConstIterator i = reg_spec.begin();
  const std::string reg_name = reg_spec.name(i);     
  const unsigned int reg_id = 9959;                   // something arbitrary
  
  Teuchos::ParameterList reg_params = reg_spec.sublist(reg_spec.name(i));
  
  // Create a rectangular region
  Teuchos::RCP<Amanzi::AmanziGeometry::Region> reg = 
    Amanzi::AmanziGeometry::createRegion(reg_spec.name(i), reg_id,
					 reg_params, &ecomm);
  
  // See if we retrieved the name and id correctly
  CHECK_EQUAL(reg->name(),reg_name);
  CHECK_EQUAL(reg->id(),reg_id);
  CHECK_EQUAL(reg_spec.isSublist(reg_spec.name(i)),true);  
  
  // Get the min-max bounds of the region from the XML specification
  Teuchos::Array<double> in_min_xyz, in_max_xyz;

  Teuchos::ParameterList::ConstIterator j = reg_params.begin();
  Teuchos::ParameterList box_params = reg_params.sublist(reg_params.name(j));
  in_min_xyz = box_params.get< Teuchos::Array<double> >("Low Coordinate");
  in_max_xyz = box_params.get< Teuchos::Array<double> >("High Coordinate");
 
  // Make sure that the region type is a BOX
  CHECK_EQUAL(reg->type(),Amanzi::AmanziGeometry::BOX);

  // Make sure that the region dimension is 3
  CHECK_EQUAL(reg->topological_dimension(),3);
  
  // See if the min-max of the region were correctly retrieved
  Amanzi::AmanziGeometry::Point pmin, pmax;
  Teuchos::RCP<const Amanzi::AmanziGeometry::RegionBox> rect =
    Teuchos::rcp_dynamic_cast<const Amanzi::AmanziGeometry::RegionBox>(reg);

  rect->corners(&pmin,&pmax);

  // Make sure we got back 3D points
  CHECK_EQUAL(pmin.dim(),3);
  CHECK_EQUAL(pmax.dim(),3);

  // Compare coordinates read from XML file and retrieved from region
  CHECK_EQUAL(pmin.x(),in_min_xyz[0]);
  CHECK_EQUAL(pmin.y(),in_min_xyz[1]);
  CHECK_EQUAL(pmin.z(),in_min_xyz[2]);
  CHECK_EQUAL(pmax.x(),in_max_xyz[0]);
  CHECK_EQUAL(pmax.y(),in_max_xyz[1]);
  CHECK_EQUAL(pmax.z(),in_max_xyz[2]);
 
  // test the functionality of the region
  std::vector<Amanzi::AmanziGeometry::Point> pin;
  pin.push_back(Amanzi::AmanziGeometry::Point(2.,3.,5.));
  pin.push_back(Amanzi::AmanziGeometry::Point(4,5,8));
  pin.push_back(Amanzi::AmanziGeometry::Point(2,3,8));
  pin.push_back(Amanzi::AmanziGeometry::Point(4,5,5));
  pin.push_back(Amanzi::AmanziGeometry::Point(2,5,5));
  pin.push_back(Amanzi::AmanziGeometry::Point(3,4,6));

  for (std::vector<Amanzi::AmanziGeometry::Point>::iterator p=pin.begin();
       p!=pin.end(); ++p) {
    CHECK(reg->inside(*p));
  }

  std::vector<Amanzi::AmanziGeometry::Point> pout;
  pin.push_back(Amanzi::AmanziGeometry::Point(3.,4.,4.9));
  pin.push_back(Amanzi::AmanziGeometry::Point(3.,4.,8.001));
  pin.push_back(Amanzi::AmanziGeometry::Point(-3,-4,-6));

  for (std::vector<Amanzi::AmanziGeometry::Point>::iterator p=pout.begin();
       p!=pout.end(); ++p) {
    CHECK(!reg->inside(*p));
  }
}  

