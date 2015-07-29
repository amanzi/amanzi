
//
// Unit test to check if a polygon region can be constructed correctly
// Author: Rao Garimella
//

#include <UnitTest++.h>

#include <iostream>


#include "../Region.hh"
#include "../PolygonRegion.hh"
#include "../RegionFactory.hh"

#include "Epetra_MpiComm.h"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Array.hpp"

#include "mpi.h"


TEST(POLYGON_REGION2)
{
  Epetra_MpiComm ecomm(MPI_COMM_WORLD);

  // read the parameter list from input file

  std::string infilename = "test/polygonregion_2D.xml";
  Teuchos::ParameterXMLFileReader xmlreader(infilename);

  Teuchos::ParameterList reg_spec(xmlreader.getParameters());


  Teuchos::ParameterList::ConstIterator i = reg_spec.begin();
  const std::string reg_name = reg_spec.name(i);     
  const unsigned int reg_id = 9959;                   // something arbitrary

  Teuchos::ParameterList reg_params = reg_spec.sublist(reg_name);

  // Create a rectangular region
  
  Amanzi::AmanziGeometry::RegionPtr reg = 
    Amanzi::AmanziGeometry::RegionFactory(reg_spec.name(i), reg_id, reg_params,
                                          2, &ecomm);
  
  // See if we retrieved the name and id correctly
  
  CHECK_EQUAL(reg->name(),reg_name);
  CHECK_EQUAL(reg->id(),reg_id);
  
  
  // Get the min-max bounds of the region from the XML specification
  
  int numpoints;
  Teuchos::Array<double> in_xyz;

  CHECK_EQUAL(reg_spec.isSublist(reg_spec.name(i)),true);

  Teuchos::ParameterList::ConstIterator j = reg_params.begin();
  Teuchos::ParameterList poly_params = reg_params.sublist(reg_params.name(j));
  numpoints = poly_params.get<int>("Number of points");
  in_xyz = poly_params.get< Teuchos::Array<double> >("Points");

  double tolerance = 1.0e-08;
  if (poly_params.isSublist("Expert Parameters")) {
    Teuchos::ParameterList expert_params = poly_params.sublist("Expert Parameters");
    tolerance = expert_params.get<double>("Tolerance");
  }
  
 
  // Make sure that the region type is a Plane

  CHECK_EQUAL(reg->type(),Amanzi::AmanziGeometry::POLYGON);
  
  // See if the parameters of the region were correctly retrieved
  
  Amanzi::AmanziGeometry::PolygonRegionPtr poly =
    dynamic_cast<Amanzi::AmanziGeometry::PolygonRegionPtr> (reg);

  int np = poly->num_points();
  CHECK_EQUAL(numpoints,np);
 
  std::vector<Amanzi::AmanziGeometry::Point> points = poly->points();
  int dim = points[0].dim();

  for (int i = 0; i < np; i++)
    for (int j = 0; j < dim; j++)
      CHECK_EQUAL(points[i][j],in_xyz[dim*i+j]);

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
  Epetra_MpiComm ecomm(MPI_COMM_WORLD);

  // read the parameter list from input file

  std::string infilename = "test/polygonregion_3D.xml";
  Teuchos::ParameterXMLFileReader xmlreader(infilename);

  Teuchos::ParameterList reg_spec(xmlreader.getParameters());


  Teuchos::ParameterList::ConstIterator i = reg_spec.begin();
  const std::string reg_name = reg_spec.name(i);     
  const unsigned int reg_id = 9959;                   // something arbitrary

  Teuchos::ParameterList reg_params = reg_spec.sublist(reg_name);

  // Create a rectangular region
  
  Amanzi::AmanziGeometry::RegionPtr reg = 
    Amanzi::AmanziGeometry::RegionFactory(reg_spec.name(i), reg_id, reg_params,
                                          3, &ecomm);
  
  // See if we retrieved the name and id correctly
  
  CHECK_EQUAL(reg->name(),reg_name);
  CHECK_EQUAL(reg->id(),reg_id);
  
  
  // Get the min-max bounds of the region from the XML specification
  
  int numpoints;
  Teuchos::Array<double> in_xyz;

  CHECK_EQUAL(reg_spec.isSublist(reg_spec.name(i)),true);

  Teuchos::ParameterList::ConstIterator j = reg_params.begin();
  Teuchos::ParameterList poly_params = reg_params.sublist(reg_params.name(j));
  numpoints = poly_params.get<int>("Number of points");
  in_xyz = poly_params.get< Teuchos::Array<double> >("Points");

  double tolerance=1.0e-08;
  if (poly_params.isSublist("Expert Parameters")) {
    Teuchos::ParameterList expert_parameters = poly_params.sublist("Expert Parameters");
    tolerance = expert_parameters.get<double>("Tolerance");
  }
  
 
  // Make sure that the region type is a Plane

  CHECK_EQUAL(reg->type(),Amanzi::AmanziGeometry::POLYGON);
  
  // See if the parameters of the region were correctly retrieved
  
  Amanzi::AmanziGeometry::PolygonRegionPtr poly =
    dynamic_cast<Amanzi::AmanziGeometry::PolygonRegionPtr> (reg);

  int np = poly->num_points();
  CHECK_EQUAL(numpoints,np);
 
  std::vector<Amanzi::AmanziGeometry::Point> points = poly->points();
  int dim = points[0].dim();

  for (int i = 0; i < np; i++)
    for (int j = 0; j < dim; j++)
      CHECK_EQUAL(points[i][j],in_xyz[dim*i+j]);

  // See if the derived parameters are sane

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
  // set to something greater than 1.0e-4 in the XML file
  testp.set(0.1,0.1,0.1001);
  CHECK(poly->inside(testp));

}  



