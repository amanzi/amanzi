/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Rao Garimella
*/

//
// Unit test to check if a plane region can be constructed correctly

#include <UnitTest++.h>

#include <iostream>


#include "AmanziComm.hh"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Array.hpp"

#include "../Point.hh"
#include "../Region.hh"
#include "../RegionPlane.hh"
#include "../RegionFactory.hh"

#include "mpi.h"


TEST(PLANE_REGION)
{
  auto ecomm = Amanzi::getDefaultComm();

  // read the parameter list from input file

  std::string infilename = "test/planeregion.xml";
  Teuchos::ParameterXMLFileReader xmlreader(infilename);

  Teuchos::ParameterList reg_spec(xmlreader.getParameters());


  Teuchos::ParameterList::ConstIterator i = reg_spec.begin();
  const std::string reg_name = reg_spec.name(i);
  const unsigned int reg_id = 9959; // something arbitrary

  Teuchos::ParameterList reg_params = reg_spec.sublist(reg_name);

  // Create a rectangular region

  Teuchos::RCP<const Amanzi::AmanziGeometry::Region> reg =
    Amanzi::AmanziGeometry::createRegion(reg_spec.name(i), reg_id, reg_params, *ecomm);

  // See if we retrieved the name and id correctly

  CHECK_EQUAL(reg->get_name(), reg_name);
  CHECK_EQUAL(reg->get_id(), reg_id);


  // Get the min-max bounds of the region from the XML specification

  Teuchos::Array<double> in_xyz, in_nrm;

  CHECK_EQUAL(reg_spec.isSublist(reg_spec.name(i)), true);

  Teuchos::ParameterList::ConstIterator j = reg_params.begin();
  Teuchos::ParameterList plane_params = reg_params.sublist(reg_params.name(j));
  in_xyz = plane_params.get<Teuchos::Array<double>>("point");
  in_nrm = plane_params.get<Teuchos::Array<double>>("normal");

  double tolerance = 1.0e-08;
  // if (plane_params.isSublist("expert parameters")) {
  //   Teuchos::ParameterList expert_params = plane_params.sublist("expert parameters");
  //   tolerance = expert_params.get<double>("tolerance");
  // }

  // Make sure that the region type is a Plane

  CHECK_EQUAL(reg->get_type(), Amanzi::AmanziGeometry::RegionType::PLANE);

  // See if the min-max of the region were correctly retrieved

  Amanzi::AmanziGeometry::Point p, n;
  Teuchos::RCP<const Amanzi::AmanziGeometry::RegionPlane> plane =
    Teuchos::rcp_dynamic_cast<const Amanzi::AmanziGeometry::RegionPlane>(reg);

  p = plane->point();
  n = plane->normal();

  double len = sqrt(in_nrm[0] * in_nrm[0] + in_nrm[1] * in_nrm[1] + in_nrm[2] * in_nrm[2]);
  CHECK_EQUAL(p.x(), in_xyz[0]);
  CHECK_EQUAL(p.y(), in_xyz[1]);
  CHECK_EQUAL(p.z(), in_xyz[2]);
  CHECK_EQUAL(n.x(), in_nrm[0] / len);
  CHECK_EQUAL(n.y(), in_nrm[1] / len);
  CHECK_EQUAL(n.z(), in_nrm[2] / len);

  int dim = p.dim();
  Amanzi::AmanziGeometry::Point p2(dim), testp(dim);

  // See if a point we know to be on the plane is considered to be "inside"

  p2.set(p.x() + 0.1, p.y() + 0.2, p.z() + 0.3); // create a point p2 that will probably
                                                 // be off the plane
  Amanzi::AmanziGeometry::Point ppvec(dim);
  ppvec = p - p2;
  double dp = ppvec * n;
  testp = p2 + dp * n;
  CHECK(plane->inside(testp));

  // Perturb the point on the plane by a small amount (smaller than
  // the comparison tolerance) and see if it is considered to be inside

  testp = p2 + (1 - (tolerance / 5)) * dp * n;
  CHECK(plane->inside(testp));

  // Perturb the point by a large amount and see if its considered to be outside

  testp = p2 + 1.1 * dp * n;
  CHECK(!plane->inside(testp));
}
