
//
// Unit test to check if a labeled set region can be constructed correctly
// Author: Rao Garimella
//

#include <UnitTest++.h>

#include <iostream>


#include "../Region.hh"
#include "../LabeledSetRegion.hh"
#include "../BoxRegion.hh"
#include "../PlaneRegion.hh"
#include "../RegionFactory.hh"
#include "../GeometricModel.hh"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Array.hpp"

#include "mpi.h"


TEST(GEOMETRIC_MODEL)
{


  // read the parameter list from input file

  std::string infilename = "manyregions.xml";
  Teuchos::ParameterXMLFileReader xmlreader(infilename);

  Teuchos::ParameterList reg_spec(xmlreader.getParameters());

  Amanzi::AmanziGeometry::GeometricModelPtr gm = new Amanzi::AmanziGeometry::GeometricModel((unsigned int)3,reg_spec);



  for (Teuchos::ParameterList::ConstIterator i = reg_spec.begin(); 
       i != reg_spec.end(); i++) {

    const std::string reg_name = reg_spec.name(i);     

    CHECK_EQUAL(reg_spec.isSublist(reg_name),true);


    Teuchos::ParameterList reg_params = reg_spec.sublist(reg_name);

    // See if the geometric model has a region by this name
  
    Amanzi::AmanziGeometry::RegionPtr reg = gm->FindRegion(reg_name);

    CHECK(reg != NULL);

    // Do their names match ?

    CHECK_EQUAL(reg->name(),reg_name);


    // Get the region info directly from the XML and compare
  
    Teuchos::ParameterList::ConstIterator j = reg_params.begin(); 

    std::string shape = reg_params.name(j);

    if (shape == "Region: Plane") {

      // Make sure that the region type is a Plane

      CHECK_EQUAL(reg->type(),Amanzi::AmanziGeometry::PLANE);
  

      // See if the point and normal of the region were correctly retrieved
  
      Teuchos::Array<double> in_xyz, in_nrm;

      Teuchos::ParameterList plane_params = reg_params.sublist(reg_params.name(j));
      in_xyz = plane_params.get< Teuchos::Array<double> >("Location");
      in_nrm = plane_params.get< Teuchos::Array<double> >("Direction");  
 
      Amanzi::AmanziGeometry::Point p, n;
      Amanzi::AmanziGeometry::PlaneRegionPtr plane =
	dynamic_cast<Amanzi::AmanziGeometry::PlaneRegionPtr> (reg);

      p = plane->point();
      n = plane->normal();
      
      CHECK_EQUAL(p.x(),in_xyz[0]);
      CHECK_EQUAL(p.y(),in_xyz[1]);
      if (p.dim() == 3)
	CHECK_EQUAL(p.z(),in_xyz[2]);
      CHECK_EQUAL(n.x(),in_nrm[0]);
      CHECK_EQUAL(n.y(),in_nrm[1]);
      if (p.dim() == 3)
	CHECK_EQUAL(n.z(),in_nrm[2]);
      
    }
    else if (shape == "Region: Box") {

      // Make sure that the region type is a BOX

      CHECK_EQUAL(reg->type(),Amanzi::AmanziGeometry::BOX);

      // Get the min-max bounds of the region from the XML specification
  
      Teuchos::Array<double> in_min_xyz, in_max_xyz;

      Teuchos::ParameterList box_params = reg_params.sublist(reg_params.name(j));
      in_min_xyz = box_params.get< Teuchos::Array<double> >("Low Coordinate");
      in_max_xyz = box_params.get< Teuchos::Array<double> >("High Coordinate");
 
      // Make sure that the region dimension is equal to the topological
      // dimension of the box
      
      CHECK_EQUAL(reg->dimension(),in_min_xyz.size());
      
      // See if the min-max of the region were correctly retrieved
      
      Amanzi::AmanziGeometry::Point pmin, pmax;
      Amanzi::AmanziGeometry::BoxRegionPtr rect =
	dynamic_cast<Amanzi::AmanziGeometry::BoxRegionPtr> (reg);
      
      rect->corners(&pmin,&pmax);
      
      // Make sure we got back 3D points
      
      CHECK_EQUAL(pmin.dim(),in_min_xyz.size());
      CHECK_EQUAL(pmax.dim(),in_min_xyz.size());
      
      // Compare coordinates read from XML file and retrieved from region
            
      CHECK_EQUAL(pmin.x(),in_min_xyz[0]);
      CHECK_EQUAL(pmin.y(),in_min_xyz[1]);      
      if (pmin.dim() == 3) {
	CHECK_EQUAL(pmin.z(),in_min_xyz[2]);
      }
      CHECK_EQUAL(pmax.x(),in_max_xyz[0]);
      CHECK_EQUAL(pmax.y(),in_max_xyz[1]);
      if (pmax.dim() == 3) {
	CHECK_EQUAL(pmax.z(),in_max_xyz[2]);
      }
    }
    else if (shape == "Region: Labeled Set") {

      Teuchos::ParameterList labset_params = reg_params.sublist(reg_params.name(j));
      std::string in_entity_str = labset_params.get< std::string >("Entity");

      // Make sure that the region type is a Labeled Set

      CHECK_EQUAL(reg->type(),Amanzi::AmanziGeometry::LABELEDSET);
  
      Amanzi::AmanziGeometry::LabeledSetRegionPtr lsreg =
	dynamic_cast<Amanzi::AmanziGeometry::LabeledSetRegionPtr> (reg);

      // Did we get the entity string right?

      CHECK_EQUAL(in_entity_str,lsreg->entity_str());
    }
  }
}  



