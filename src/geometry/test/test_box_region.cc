#include <UnitTest++.h>

#include <iostream>


#include "../Region.hh"
#include "../BoxRegion.hh"
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
    
  std::cout << reg_params;

  // Create a rectangular region
  
  Amanzi::AmanziGeometry::RegionPtr reg = 
    Amanzi::AmanziGeometry::RegionFactory(reg_name, reg_id, reg_params, &ecomm);
  
  // See if we retrieved the name and id correctly
  
  CHECK_EQUAL(reg->name(),reg_name);
  CHECK_EQUAL(reg->id(),reg_id);
  
  
  // Get the min-max bounds of the region from the XML specification
  
  Teuchos::Array<double> in_min_xyz, in_max_xyz;

  CHECK_EQUAL(reg_spec.isSublist(reg_spec.name(i)),true);

  Teuchos::ParameterList::ConstIterator j = reg_params.begin();
  Teuchos::ParameterList box_params = reg_params.sublist(reg_params.name(j));
  in_min_xyz = box_params.get< Teuchos::Array<double> >("Low Coordinate");
  in_max_xyz = box_params.get< Teuchos::Array<double> >("High Coordinate");

  
 
  // Make sure that the region type is a BOX

  CHECK_EQUAL(reg->type(),Amanzi::AmanziGeometry::BOX);

  // Make sure that the region dimension is 2

  CHECK_EQUAL(reg->dimension(),2);
  
  // See if the min-max of the region were correctly retrieved
  
  Amanzi::AmanziGeometry::Point pmin, pmax;
  Amanzi::AmanziGeometry::BoxRegionPtr rect =
    dynamic_cast<Amanzi::AmanziGeometry::BoxRegionPtr> (reg);

  rect->corners(&pmin,&pmax);

  // Make sure we got back 2D points

  CHECK_EQUAL(pmin.dim(),2);
  CHECK_EQUAL(pmax.dim(),2);

  // Compare coordinates read from XML file and retrieved from region

  CHECK_EQUAL(pmin.x(),in_min_xyz[0]);
  CHECK_EQUAL(pmin.y(),in_min_xyz[1]);
  CHECK_EQUAL(pmax.x(),in_max_xyz[0]);
  CHECK_EQUAL(pmax.y(),in_max_xyz[1]);
 
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
  
  Amanzi::AmanziGeometry::RegionPtr reg = 
    Amanzi::AmanziGeometry::RegionFactory(reg_spec.name(i), reg_id, reg_params, &ecomm);
  
  // See if we retrieved the name and id correctly
  
  CHECK_EQUAL(reg->name(),reg_name);
  CHECK_EQUAL(reg->id(),reg_id);
  
  
  // Get the min-max bounds of the region from the XML specification
  
  Teuchos::Array<double> in_min_xyz, in_max_xyz;

  CHECK_EQUAL(reg_spec.isSublist(reg_spec.name(i)),true);

  Teuchos::ParameterList::ConstIterator j = reg_params.begin();
  Teuchos::ParameterList box_params = reg_params.sublist(reg_params.name(j));
  in_min_xyz = box_params.get< Teuchos::Array<double> >("Low Coordinate");
  in_max_xyz = box_params.get< Teuchos::Array<double> >("High Coordinate");

  
 
  // Make sure that the region type is a BOX

  CHECK_EQUAL(reg->type(),Amanzi::AmanziGeometry::BOX);

  // Make sure that the region dimension is 3

  CHECK_EQUAL(reg->dimension(),3);
  
  // See if the min-max of the region were correctly retrieved
  
  Amanzi::AmanziGeometry::Point pmin, pmax;
  Amanzi::AmanziGeometry::BoxRegionPtr rect =
    dynamic_cast<Amanzi::AmanziGeometry::BoxRegionPtr> (reg);

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
 
}  

