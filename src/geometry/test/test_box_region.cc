#include <UnitTest++.h>

#include <iostream>
#include <vector>

#include "../Region.hh"
#include "../RegionBox.hh"
#include "../RegionBoxVolumeFractions.hh"
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
  in_min_xyz = box_params.get< Teuchos::Array<double> >("low coordinate");
  in_max_xyz = box_params.get< Teuchos::Array<double> >("high coordinate");
 
  // Make sure that the region type is a BOX
  CHECK_EQUAL(reg->type(),Amanzi::AmanziGeometry::BOX);

  // Make sure that the region dimension is 2
  CHECK_EQUAL(reg->manifold_dimension(),2);
  
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
  in_min_xyz = box_params.get< Teuchos::Array<double> >("low coordinate");
  in_max_xyz = box_params.get< Teuchos::Array<double> >("high coordinate");
 
  // Make sure that the region type is a BOX
  CHECK_EQUAL(reg->type(),Amanzi::AmanziGeometry::BOX);

  // Make sure that the region dimension is 3
  CHECK_EQUAL(reg->manifold_dimension(),3);
  
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


TEST(BOXREGION_VOFS_2D_INTERSECTION)
{
  Amanzi::AmanziGeometry::Point v1(2), v2(2), v3(2), v4(2), vv(2);
  std::vector<Amanzi::AmanziGeometry::Point> xy1, xy2, xy3;

  v1.set(0.0, 0.0);
  v2.set(1.0, 0.0);
  v3.set(0.0, 1.0);
  v4.set(1.0, 1.0);

  xy2.push_back(v1);
  xy2.push_back(v2);
  xy2.push_back(v4);
  xy2.push_back(v3);

  int n(0), sizes[6] = {3, 5, 5, 4, 4, 0};
  for (double d = 0.0; d <= 1.0; d += 0.2) {
    vv.set(d, d);
    xy1.clear();
    xy1.push_back(vv + v1);
    xy1.push_back(vv + v2);
    xy1.push_back(vv + v3);
    std::cout << "\nshift: " << xy1[0] << std::endl;

    Amanzi::AmanziGeometry::IntersectConvexPolygons(xy1, xy2, xy3);

    for (int i = 0; i < xy3.size(); ++i) {
      std::cout << i << " xy=" << xy3[i] << std::endl;
    }
    CHECK(xy3.size() == sizes[n++]);
  }
}


TEST(BOXREGION_VOFS_2D_AREA)
{
  Epetra_MpiComm ecomm(MPI_COMM_WORLD);

  // read the parameter list from input file
  std::string infilename = "test/boxregion_vofs.xml";
  Teuchos::ParameterXMLFileReader xmlreader(infilename);
  Teuchos::ParameterList reg_spec(xmlreader.getParameters());

  // create a rectangular region
  Teuchos::ParameterList::ConstIterator i = reg_spec.begin();
  std::string reg_name = reg_spec.name(i);     
  unsigned int reg_id = 9959;  // something arbitrary
  Teuchos::ParameterList reg_params = reg_spec.sublist(reg_name);
    
  Teuchos::RCP<Amanzi::AmanziGeometry::Region> reg = 
    Amanzi::AmanziGeometry::createRegion(reg_name, reg_id, reg_params, &ecomm);
  
  Amanzi::AmanziGeometry::Point v1(2), v2(2), v3(2), vv(2);
  std::vector<Amanzi::AmanziGeometry::Point> polygon;

  v1.set(0.0, 0.0);
  v2.set(1.0, 0.0);
  v3.set(0.0, 1.0);

  int n(0);
  double area_exact[5] = {0.5, 0.46, 0.34, 0.16, 0.04};
  for (double d = 0.0; d <= 0.8; d += 0.2) {
    vv.set(d, d);
    polygon.clear();
    polygon.push_back(vv + v1);
    polygon.push_back(vv + v2);
    polygon.push_back(vv + v3);

    double area = reg->intersect(polygon);
    CHECK_CLOSE(area, area_exact[n++], 1e-6);
  }
}


TEST(BOXREGION_VOFS_3D_INTERSECTION)
{
  using namespace Amanzi::AmanziGeometry;

  Point v1(3);
  std::vector<Point> xy1, xy3;
  std::vector<std::vector<int> > faces1(4), faces3;
  std::vector<std::pair<Point, Point> > xy2;

  xy2.push_back(std::make_pair(Point(0.0, 0.0, 0.0), Point(-1.0, 0.0, 0.0)));
  xy2.push_back(std::make_pair(Point(0.0, 0.0, 0.0), Point(0.0, -1.0, 0.0)));
  xy2.push_back(std::make_pair(Point(0.0, 0.0, 0.0), Point(0.0, 0.0, -1.0)));

  xy2.push_back(std::make_pair(Point(1.0, 1.0, 1.0), Point(1.0, 0.0, 0.0)));
  xy2.push_back(std::make_pair(Point(1.0, 1.0, 1.0), Point(0.0, 1.0, 0.0)));
  xy2.push_back(std::make_pair(Point(1.0, 1.0, 1.0), Point(0.0, 0.0, 1.0)));

  xy1.push_back(Point(0.0, 0.0, 0.0));
  xy1.push_back(Point(1.0, 0.0, 0.0));
  xy1.push_back(Point(0.0, 1.0, 0.0));
  xy1.push_back(Point(0.0, 0.0, 1.0));

  faces1[0].push_back(0);
  faces1[0].push_back(2);
  faces1[0].push_back(1);

  faces1[1].push_back(0);
  faces1[1].push_back(1);
  faces1[1].push_back(3);

  faces1[2].push_back(0);
  faces1[2].push_back(3);
  faces1[2].push_back(2);

  faces1[3].push_back(1);
  faces1[3].push_back(2);
  faces1[3].push_back(3);

  Amanzi::AmanziGeometry::IntersectConvexPolyhedra(xy1, faces1, xy2, xy3, faces3);
}

