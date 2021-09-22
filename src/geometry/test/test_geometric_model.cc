
//
// Unit test to check if a labeled set region can be constructed correctly
// Author: Rao Garimella
//

#include <UnitTest++.h>

#include <iostream>

#include "AmanziComm.hh"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Array.hpp"

#include "Region.hh"
#include "RegionLabeledSet.hh"
#include "RegionBox.hh"
#include "RegionPlane.hh"
#include "RegionFactory.hh"
#include "GeometricModel.hh"

TEST(GEOMETRIC_MODEL)
{

  auto ecomm = Amanzi::getDefaultComm();

  // read the parameter list from input file

  std::string infilename = "test/manyregions.xml";
  Teuchos::ParameterXMLFileReader xmlreader(infilename);

  Teuchos::ParameterList reg_spec(xmlreader.getParameters());

  Amanzi::AmanziGeometry::GeometricModel gm((unsigned int)3, reg_spec, *ecomm);

  for (Teuchos::ParameterList::ConstIterator i = reg_spec.begin();
       i != reg_spec.end(); i++) {
    const std::string reg_name = reg_spec.name(i);
    CHECK_EQUAL(reg_spec.isSublist(reg_name),true);

    Teuchos::ParameterList reg_params = reg_spec.sublist(reg_name);

    // See if the geometric model has a region by this name
    Teuchos::RCP<const Amanzi::AmanziGeometry::Region> reg = gm.FindRegion(reg_name);

    CHECK(reg.get());

    // Do their names match ?
    CHECK_EQUAL(reg->get_name(),reg_name);

    // Get the region info directly from the XML and compare
    Teuchos::ParameterList::ConstIterator j = reg_params.begin();

    std::string shape = reg_params.name(j);

    if (shape == "region: plane") {
      // Make sure that the region type is a Plane
      CHECK_EQUAL( Amanzi::AmanziGeometry::RegionType::PLANE,  reg->get_type());

      // See if the point and normal of the region were correctly retrieved
      Teuchos::Array<double> in_xyz, in_nrm;

      Teuchos::ParameterList plane_params = reg_params.sublist(reg_params.name(j));
      in_xyz = plane_params.get< Teuchos::Array<double> >("point");
      in_nrm = plane_params.get< Teuchos::Array<double> >("normal");

      Amanzi::AmanziGeometry::Point p, n;
      Teuchos::RCP<const Amanzi::AmanziGeometry::RegionPlane> plane =
        Teuchos::rcp_dynamic_cast<const Amanzi::AmanziGeometry::RegionPlane>(reg);

      p = plane->point();
      n = plane->normal();

      CHECK_EQUAL(p.x(),in_xyz[0]);
      CHECK_EQUAL(p.y(),in_xyz[1]);
      if (p.dim() == 3)
        CHECK_EQUAL(p.z(),in_xyz[2]);
      double len = in_nrm[0]*in_nrm[0]+in_nrm[1]*in_nrm[1];
      if (p.dim() == 3) len += in_nrm[2]*in_nrm[2];
      len = sqrt(len);
      CHECK_EQUAL(n.x(),in_nrm[0]/len);
      CHECK_EQUAL(n.y(),in_nrm[1]/len);
      if (p.dim() == 3)
        CHECK_EQUAL(n.z(),in_nrm[2]/len);

    } else if (shape == "region: box") {
      // Make sure that the region type is a BOX
      CHECK_EQUAL( reg->get_type(), Amanzi::AmanziGeometry::RegionType::BOX);

      // Get the min-max bounds of the region from the XML specification
      Teuchos::Array<double> in_min_xyz, in_max_xyz;

      Teuchos::ParameterList box_params = reg_params.sublist(reg_params.name(j));
      in_min_xyz = box_params.get< Teuchos::Array<double> >("low coordinate");
      in_max_xyz = box_params.get< Teuchos::Array<double> >("high coordinate");

      // Make sure that the region dimension is equal to the topological
      // dimension of the box
      CHECK_EQUAL(reg->get_manifold_dimension(),in_min_xyz.size());

      // See if the min-max of the region were correctly retrieved
      Amanzi::AmanziGeometry::Point pmin, pmax;
      auto rect = Teuchos::rcp_dynamic_cast<const Amanzi::AmanziGeometry::RegionBox>(reg);

      pmin = rect->point0();
      pmax = rect->point1();

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

    } else if (shape == "region: labeled set") {
      Teuchos::ParameterList labset_params = reg_params.sublist(reg_params.name(j));
      std::string in_entity_str = labset_params.get< std::string >("entity");

      // Make sure that the region type is a Labeled Set
      CHECK_EQUAL(reg->get_type(), Amanzi::AmanziGeometry::RegionType::LABELEDSET);

      Teuchos::RCP<const Amanzi::AmanziGeometry::RegionLabeledSet> lsreg =
        Teuchos::rcp_dynamic_cast<const Amanzi::AmanziGeometry::RegionLabeledSet>(reg);

      // Did we get the entity string right?
      CHECK_EQUAL(in_entity_str,lsreg->entity_str());
    }
  }
}
