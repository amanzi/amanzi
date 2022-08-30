//
// Unit test to check if a polygon region can be constructed correctly
// Author: Rao Garimella
//

#include <iostream>

#include "mpi.h"
#include "UnitTest++.h"
#include "AmanziComm.hh"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Array.hpp"

#include "../Point.hh"
#include "../Region.hh"
#include "../RegionPolygon.hh"
#include "../RegionFactory.hh"


TEST(POLYGON_REGION2)
{
  auto ecomm = Amanzi::getDefaultComm();

  // read the parameter list from input file
  std::string infilename = "test/polygonregion_2D.xml";
  Teuchos::ParameterXMLFileReader xmlreader(infilename);

  Teuchos::ParameterList reg_spec(xmlreader.getParameters());


  Teuchos::ParameterList::ConstIterator i = reg_spec.begin();
  const std::string reg_name = reg_spec.name(i);     
  const unsigned int reg_id = 9959;                   // something arbitrary

  Teuchos::ParameterList reg_params = reg_spec.sublist(reg_name);

  // Create a rectangular region
  Teuchos::RCP<const Amanzi::AmanziGeometry::Region> reg = 
    Amanzi::AmanziGeometry::createRegion(reg_spec.name(i), reg_id,
					 reg_params, *ecomm);
  
  // See if we retrieved the name and id correctly
  CHECK_EQUAL(reg->get_name(),reg_name);
  CHECK_EQUAL(reg->get_id(),reg_id);
  
  // Get the min-max bounds of the region from the XML specification
  int numpoints;
  Teuchos::Array<double> in_xyz;
  CHECK_EQUAL(reg_spec.isSublist(reg_spec.name(i)),true);

  Teuchos::ParameterList poly_params = reg_params.sublist(reg_params.name(reg_params.begin()));
  numpoints = poly_params.get<int>("number of points");
  in_xyz = poly_params.get< Teuchos::Array<double> >("points");

  double tolerance = 1.0e-08;
  // if (poly_params.isSublist("expert parameters")) {
  //   Teuchos::ParameterList expert_params = poly_params.sublist("expert parameters");
  //   tolerance = expert_params.get<double>("tolerance");
  // }
  
  // Make sure that the region type is a Plane
  CHECK_EQUAL(reg->get_type(),Amanzi::AmanziGeometry::RegionType::POLYGON);
  
  // See if the parameters of the region were correctly retrieved
  Teuchos::RCP<const Amanzi::AmanziGeometry::RegionPolygon> poly =
    Teuchos::rcp_dynamic_cast<const Amanzi::AmanziGeometry::RegionPolygon>(reg);

  int np = poly->size();
  CHECK_EQUAL(numpoints,np);

  int lcv=0;
  int dim = 0;
  for (Amanzi::AmanziGeometry::RegionPolygon::PointIterator p=poly->begin();
       p!=poly->end(); ++lcv,++p) {
    dim = p->dim();
    for (int j = 0; j < dim; j++)
      CHECK_EQUAL((*p)[j],in_xyz[dim*lcv+j]);
  }

  // See if the derived parameters are sane

  Amanzi::AmanziGeometry::Point normal = poly->normal();
  CHECK_CLOSE(normal[0],sqrt(0.5),1.0e-06);
  CHECK_CLOSE(normal[1],sqrt(0.5),1.0e-06);


  // See if a point we know is considered to be inside

  Amanzi::AmanziGeometry::Point testp(dim);
  testp.set(0.0,0.0);
  CHECK(poly->inside(testp));

  // Check a point that is on the boundary of the polygon
  testp.set(0.5,-0.5);
  CHECK(poly->inside(testp));
    
  // Check a point we know to be off the plane
  testp.set(0.0,0.1);
  CHECK(!poly->inside(testp));

  // Check a point we know to be on the plane but outside the polygon
  testp.set(0.9,0.9);
  CHECK(!poly->inside(testp));

  // Check a point we know to be very close to the plane (within tolerance)
  testp.set(0.0,tolerance/2.0);
  CHECK(poly->inside(testp));
}  

TEST(POLYGON_REGION3)
{
  auto ecomm = Amanzi::getDefaultComm();

  // read the parameter list from input file

  std::string infilename = "test/polygonregion_3D.xml";
  Teuchos::ParameterXMLFileReader xmlreader(infilename);

  Teuchos::ParameterList reg_spec(xmlreader.getParameters());


  Teuchos::ParameterList::ConstIterator i = reg_spec.begin();
  const std::string reg_name = reg_spec.name(i);     
  const unsigned int reg_id = 9959;                   // something arbitrary

  Teuchos::ParameterList reg_params = reg_spec.sublist(reg_name);

  // Create a rectangular region
  
  Teuchos::RCP<const Amanzi::AmanziGeometry::Region> reg = 
    Amanzi::AmanziGeometry::createRegion(reg_spec.name(i), reg_id,
					 reg_params, *ecomm);
  
  // See if we retrieved the name and id correctly
  CHECK_EQUAL(reg->get_name(),reg_name);
  CHECK_EQUAL(reg->get_id(),reg_id);
  CHECK_EQUAL(reg_spec.isSublist(reg_spec.name(i)),true);
  
  // Get the min-max bounds of the region from the XML specification
  int numpoints;
  Teuchos::Array<double> in_xyz;

  Teuchos::ParameterList poly_params = reg_params.sublist(reg_params.name(reg_params.begin()));
  numpoints = poly_params.get<int>("number of points");
  in_xyz = poly_params.get< Teuchos::Array<double> >("points");
 
  // Make sure that the region type is a Plane
  CHECK_EQUAL(reg->get_type(),Amanzi::AmanziGeometry::RegionType::POLYGON);
  
  // See if the parameters of the region were correctly retrieved
  Teuchos::RCP<const Amanzi::AmanziGeometry::RegionPolygon> poly =
    Teuchos::rcp_dynamic_cast<const Amanzi::AmanziGeometry::RegionPolygon>(reg);

  int np = poly->size();
  CHECK_EQUAL(numpoints,np);
 
  int lcv=0;
  int dim = 0;
  for (Amanzi::AmanziGeometry::RegionPolygon::PointIterator p=poly->begin();
       p!=poly->end(); ++lcv,++p) {
    dim = p->dim();
    for (int j = 0; j < dim; j++)
      CHECK_EQUAL((*p)[j],in_xyz[dim*lcv+j]);
  }

  Amanzi::AmanziGeometry::Point normal = poly->normal();
  CHECK_CLOSE(normal[0],0.0,1.0e-06);
  CHECK_CLOSE(normal[1],-sqrt(0.5),1.0e-06);
  CHECK_CLOSE(normal[2],sqrt(0.5),1.0e-06);


  // See if a point we know is considered to be inside
  Amanzi::AmanziGeometry::Point testp(dim);
  testp.set(0.1,0.1,0.1);
  CHECK(poly->inside(testp));

  // Check a point known to be on a vertex of the polygon
  testp.set(0.5,-0.5,-0.5);     // passes
  testp.set(-0.5,-0.5,-0.5);    // fails
  CHECK(poly->inside(testp));

  // Check a point known to be on an edge of the polygon
  testp.set(0.1,-0.5,-0.5);    // passes
  testp.set(-0.5,0.1,0.1);     // fails
  CHECK(poly->inside(testp));

  // Check a point along the infinite line of an edge of the polygon
  // but outside the edge

  testp.set(-0.9,-0.5,-0.5);
  CHECK(!poly->inside(testp));

  // Check a point we know to be off the plane
  testp.set(0.1,0.1,-0.9);
  CHECK(!poly->inside(testp));

  // Check a point we know to be on the plane but outside the polygon
  testp.set(1.0,0.0,0.0);
  CHECK(!poly->inside(testp));

  // Check a point that is close to the plane - assume tolerance has been
  // set to something greater than 1.0e-8 in the XML file
  testp.set(0.1,0.1,0.10000001);
  CHECK(poly->inside(testp));

}  



