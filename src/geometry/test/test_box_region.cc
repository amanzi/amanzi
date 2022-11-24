#include <UnitTest++.h>

#include <iostream>
#include <vector>

#include "../Region.hh"
#include "../RegionBox.hh"
#include "../RegionBoxVolumeFractions.hh"
#include "../RegionFactory.hh"

#include "AmanziComm.hh"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Array.hpp"

#include "mpi.h"

TEST(BOX_REGION_2D)
{
  auto ecomm = Amanzi::getDefaultComm();

  // read the parameter list from input file
  std::string infilename = "test/boxregion_2D.xml";
  Teuchos::ParameterXMLFileReader xmlreader(infilename);

  Teuchos::ParameterList reg_spec(xmlreader.getParameters());

  Teuchos::ParameterList::ConstIterator i = reg_spec.begin();
  const std::string reg_name = reg_spec.name(i);
  const unsigned int reg_id = 9959; // something arbitrary

  Teuchos::ParameterList reg_params = reg_spec.sublist(reg_name);

  // Create a rectangular region
  auto reg = Amanzi::AmanziGeometry::createRegion(reg_name, reg_id, reg_params, *ecomm);

  // See if we retrieved the name and id correctly
  CHECK_EQUAL(reg->get_name(), reg_name);
  CHECK_EQUAL(reg->get_id(), reg_id);
  CHECK_EQUAL(reg_spec.isSublist(reg_spec.name(i)), true);

  // Get the min-max bounds of the region from the XML specification
  Teuchos::Array<double> in_min_xyz, in_max_xyz;

  Teuchos::ParameterList::ConstIterator j = reg_params.begin();
  Teuchos::ParameterList box_params = reg_params.sublist(reg_params.name(j));
  in_min_xyz = box_params.get<Teuchos::Array<double>>("low coordinate");
  in_max_xyz = box_params.get<Teuchos::Array<double>>("high coordinate");

  // Make sure that the region type is a BOX
  CHECK_EQUAL(reg->get_type(), Amanzi::AmanziGeometry::RegionType::BOX);

  // Make sure that the region dimension is 2
  CHECK_EQUAL(reg->get_manifold_dimension(), 2);

  // See if the min-max of the region were correctly retrieved
  Amanzi::AmanziGeometry::Point pmin, pmax;
  auto rect = Teuchos::rcp_dynamic_cast<const Amanzi::AmanziGeometry::RegionBox>(reg);

  // Make sure we got back 2D points
  pmin = rect->point0();
  pmax = rect->point1();
  CHECK_EQUAL(pmin.dim(), 2);
  CHECK_EQUAL(pmax.dim(), 2);

  // Compare coordinates read from XML file and retrieved from region
  CHECK_EQUAL(pmin.x(), in_min_xyz[0]);
  CHECK_EQUAL(pmin.y(), in_max_xyz[1]);
  CHECK_EQUAL(pmax.x(), in_max_xyz[0]);
  CHECK_EQUAL(pmax.y(), in_min_xyz[1]);

  // test the functionality of the region
  std::vector<Amanzi::AmanziGeometry::Point> pin;
  pin.push_back(Amanzi::AmanziGeometry::Point(9., 8.));
  pin.push_back(Amanzi::AmanziGeometry::Point(11., 8.));
  pin.push_back(Amanzi::AmanziGeometry::Point(11., 1.5));
  pin.push_back(Amanzi::AmanziGeometry::Point(9., 1.5));
  pin.push_back(Amanzi::AmanziGeometry::Point(10., 5.));

  for (std::vector<Amanzi::AmanziGeometry::Point>::iterator p = pin.begin(); p != pin.end(); ++p) {
    CHECK(reg->inside(*p));
  }

  std::vector<Amanzi::AmanziGeometry::Point> pout;
  pin.push_back(Amanzi::AmanziGeometry::Point(9.9, 8.));
  pin.push_back(Amanzi::AmanziGeometry::Point(11, 7.9));

  for (std::vector<Amanzi::AmanziGeometry::Point>::iterator p = pout.begin(); p != pout.end();
       ++p) {
    CHECK(!reg->inside(*p));
  }
}


TEST(BOX_REGION_3D)
{
  auto ecomm = Amanzi::getDefaultComm();

  // read the parameter list from input file

  std::string infilename = "test/boxregion_3D.xml";
  Teuchos::ParameterXMLFileReader xmlreader(infilename);

  Teuchos::ParameterList reg_spec(xmlreader.getParameters());


  Teuchos::ParameterList::ConstIterator i = reg_spec.begin();
  const std::string reg_name = reg_spec.name(i);
  const unsigned int reg_id = 9959; // something arbitrary

  Teuchos::ParameterList reg_params = reg_spec.sublist(reg_spec.name(i));

  // Create a rectangular region
  Teuchos::RCP<Amanzi::AmanziGeometry::Region> reg =
    Amanzi::AmanziGeometry::createRegion(reg_spec.name(i), reg_id, reg_params, *ecomm);

  // See if we retrieved the name and id correctly
  CHECK_EQUAL(reg->get_name(), reg_name);
  CHECK_EQUAL(reg->get_id(), reg_id);
  CHECK_EQUAL(reg_spec.isSublist(reg_spec.name(i)), true);

  // Get the min-max bounds of the region from the XML specification
  Teuchos::Array<double> in_min_xyz, in_max_xyz;

  Teuchos::ParameterList::ConstIterator j = reg_params.begin();
  Teuchos::ParameterList box_params = reg_params.sublist(reg_params.name(j));
  in_min_xyz = box_params.get<Teuchos::Array<double>>("low coordinate");
  in_max_xyz = box_params.get<Teuchos::Array<double>>("high coordinate");

  // Make sure that the region type is a BOX
  CHECK_EQUAL(reg->get_type(), Amanzi::AmanziGeometry::RegionType::BOX);

  // Make sure that the region dimension is 3
  CHECK_EQUAL(reg->get_manifold_dimension(), 3);

  // See if the min-max of the region were correctly retrieved
  Amanzi::AmanziGeometry::Point pmin, pmax;
  auto rect = Teuchos::rcp_dynamic_cast<const Amanzi::AmanziGeometry::RegionBox>(reg);

  pmin = rect->point0();
  pmax = rect->point1();

  // Make sure we got back 3D points
  CHECK_EQUAL(pmin.dim(), 3);
  CHECK_EQUAL(pmax.dim(), 3);

  // Compare coordinates read from XML file and retrieved from region
  CHECK_EQUAL(pmin.x(), in_min_xyz[0]);
  CHECK_EQUAL(pmin.y(), in_min_xyz[1]);
  CHECK_EQUAL(pmin.z(), in_min_xyz[2]);
  CHECK_EQUAL(pmax.x(), in_max_xyz[0]);
  CHECK_EQUAL(pmax.y(), in_max_xyz[1]);
  CHECK_EQUAL(pmax.z(), in_max_xyz[2]);

  // test the functionality of the region
  std::vector<Amanzi::AmanziGeometry::Point> pin;
  pin.push_back(Amanzi::AmanziGeometry::Point(2., 3., 5.));
  pin.push_back(Amanzi::AmanziGeometry::Point(4, 5, 8));
  pin.push_back(Amanzi::AmanziGeometry::Point(2, 3, 8));
  pin.push_back(Amanzi::AmanziGeometry::Point(4, 5, 5));
  pin.push_back(Amanzi::AmanziGeometry::Point(2, 5, 5));
  pin.push_back(Amanzi::AmanziGeometry::Point(3, 4, 6));

  for (std::vector<Amanzi::AmanziGeometry::Point>::iterator p = pin.begin(); p != pin.end(); ++p) {
    CHECK(reg->inside(*p));
  }

  std::vector<Amanzi::AmanziGeometry::Point> pout;
  pin.push_back(Amanzi::AmanziGeometry::Point(3., 4., 4.9));
  pin.push_back(Amanzi::AmanziGeometry::Point(3., 4., 8.001));
  pin.push_back(Amanzi::AmanziGeometry::Point(-3, -4, -6));

  for (std::vector<Amanzi::AmanziGeometry::Point>::iterator p = pout.begin(); p != pout.end();
       ++p) {
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

  int n(0), sizes[6] = { 3, 5, 5, 4, 4, 0 };
  for (double d = 0.0; d <= 1.0; d += 0.2) {
    vv.set(d, d);
    xy1.clear();
    xy1.push_back(vv + v1);
    xy1.push_back(vv + v2);
    xy1.push_back(vv + v3);
    std::cout << "\nshift: " << xy1[0] << std::endl;

    Amanzi::AmanziGeometry::IntersectConvexPolygons(xy1, xy2, xy3);

    for (int i = 0; i < xy3.size(); ++i) { std::cout << i << " xy=" << xy3[i] << std::endl; }
    CHECK(xy3.size() == sizes[n++]);
  }
}

TEST(BOXREGION_VOFS_2D_AREA)
{
  using namespace Amanzi::AmanziGeometry;

  auto ecomm = Amanzi::getDefaultComm();

  // read the parameter list from input file
  std::string infilename = "test/boxregion_vofs.xml";
  Teuchos::ParameterXMLFileReader xmlreader(infilename);
  Teuchos::ParameterList reg_spec(xmlreader.getParameters());

  // create a rectangular region
  Teuchos::ParameterList::ConstIterator i = reg_spec.begin();
  std::string reg_name = reg_spec.name(i);
  unsigned int reg_id = 9959; // something arbitrary
  Teuchos::ParameterList reg_params = reg_spec.sublist(reg_name);

  Teuchos::RCP<Amanzi::AmanziGeometry::Region> reg =
    Amanzi::AmanziGeometry::createRegion(reg_name, reg_id, reg_params, *ecomm);

  Amanzi::AmanziGeometry::Point v1(2), v2(2), v3(2), vv(2);
  std::vector<Amanzi::AmanziGeometry::Point> polygon;

  v1.set(0.0, 0.0);
  v2.set(1.0, 0.0);
  v3.set(0.0, 1.0);

  int n(0);
  double area_exact[5] = { 0.5, 0.46, 0.34, 0.16, 0.04 };
  for (double d = 0.0; d <= 0.8; d += 0.2) {
    vv.set(d, d);
    polygon.clear();
    polygon.push_back(vv + v1);
    polygon.push_back(vv + v2);
    polygon.push_back(vv + v3);

    double area = reg->intersect(polygon);
    CHECK_CLOSE(area_exact[n++], area, 1e-6);
  }
}


TEST(BOXREGION_VOFS_3D_INTERSECTION)
{
  using namespace Amanzi::AmanziGeometry;

  std::cout << "\nIntersection of a moving tet with the unit box\n\n";

  Point v1(3), vv(3);
  std::vector<Point> xyz1, xyz3;
  std::vector<std::vector<int>> faces1(4), faces3;
  std::vector<std::pair<Point, Point>> xyz2;

  xyz2.push_back(std::make_pair(Point(0.0, 0.0, 0.0), Point(-1.0, 0.0, 0.0)));
  xyz2.push_back(std::make_pair(Point(0.0, 0.0, 0.0), Point(0.0, -1.0, 0.0)));
  xyz2.push_back(std::make_pair(Point(0.0, 0.0, 0.0), Point(0.0, 0.0, -1.0)));

  xyz2.push_back(std::make_pair(Point(1.0, 1.0, 1.0), Point(1.0, 0.0, 0.0)));
  xyz2.push_back(std::make_pair(Point(1.0, 1.0, 1.0), Point(0.0, 1.0, 0.0)));
  xyz2.push_back(std::make_pair(Point(1.0, 1.0, 1.0), Point(0.0, 0.0, 1.0)));

  int n(0), sizes[7] = { 4, 7, 7, 7, 6, 0, 0 };
  for (double d = 0.0; d <= 1.21; d += 0.2) {
    vv.set(d, d, d);
    std::cout << "\nShift: " << vv << std::endl;
    xyz1.clear();
    xyz1.push_back(vv + Point(0.0, 0.0, 0.0));
    xyz1.push_back(vv + Point(1.0, 0.0, 0.0));
    xyz1.push_back(vv + Point(0.0, 1.0, 0.0));
    xyz1.push_back(vv + Point(0.0, 0.0, 1.0));

    for (int i = 0; i < 4; ++i) faces1[i].clear();
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

    Amanzi::AmanziGeometry::IntersectConvexPolyhedra(xyz1, faces1, xyz2, xyz3, faces3);

    int nfaces3(faces3.size());
    std::cout << "Total number of faces: " << nfaces3 << std::endl;
    for (int i = 0; i < nfaces3; ++i) {
      int nnodes(faces3[i].size());
      for (int m = 0; m < nnodes; ++m) std::cout << faces3[i][m] << " ";
      std::cout << std::endl;
    }
    std::cout << "Total number of vertices: " << xyz3.size() << std::endl;

    CHECK(nfaces3 == sizes[n++]);
  }
}


TEST(BOXREGION_VOFS_3D_VOLUME)
{
  using namespace Amanzi::AmanziGeometry;

  std::cout << "\nVolume of intersection of a moving tet with the unit box\n\n";

  auto ecomm = Amanzi::getDefaultComm();

  // read the parameter list from input file
  std::string infilename = "test/boxregion_vofs_3D.xml";
  Teuchos::ParameterXMLFileReader xmlreader(infilename);
  Teuchos::ParameterList reg_spec(xmlreader.getParameters());

  // create a rectangular region
  std::string reg_name = reg_spec.name(reg_spec.begin());
  unsigned int reg_id = 9959; // something arbitrary
  Teuchos::ParameterList reg_params = reg_spec.sublist(reg_name);

  Teuchos::RCP<Amanzi::AmanziGeometry::Region> reg =
    Amanzi::AmanziGeometry::createRegion(reg_name, reg_id, reg_params, *ecomm);

  std::vector<Point> xyz;
  std::vector<std::vector<int>> faces(4);

  // double volume_exact[5] = {0.5, 0.46, 0.34, 0.16, 0.04};
  for (double d = 0.0; d <= 0.8; d += 0.2) {
    Point vv(d, d, d);
    xyz.clear();
    xyz.push_back(vv + Point(0.0, 0.0, 0.0));
    xyz.push_back(vv + Point(1.0, 0.0, 0.0));
    xyz.push_back(vv + Point(0.0, 1.0, 0.0));
    xyz.push_back(vv + Point(0.0, 0.0, 1.0));

    for (int i = 0; i < 4; ++i) faces[i].clear();
    faces[0].push_back(0);
    faces[0].push_back(2);
    faces[0].push_back(1);

    faces[1].push_back(0);
    faces[1].push_back(1);
    faces[1].push_back(3);

    faces[2].push_back(0);
    faces[2].push_back(3);
    faces[2].push_back(2);

    faces[3].push_back(1);
    faces[3].push_back(2);
    faces[3].push_back(3);

    double volume = reg->intersect(xyz, faces);
    std::cout << "volume=" << volume << std::endl;
  }
}
