/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/*
  Nonmember function for creating regions.

  Copyright 2010-2013 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Rao Garimella
           Ethan Coon (ecoon@lanl.gov)
*/

#include <iostream>
#include <sstream>

#include "Teuchos_Array.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_StrUtils.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "dbc.hh"
#include "errors.hh"

#include "HDF5Reader.hh"

#include "Region.hh"
#include "RegionAll.hh"
#include "RegionBoundary.hh"
#include "RegionBox.hh"
#include "RegionBoxVolumeFractions.hh"
#include "RegionCylinder.hh"
#include "RegionEnumerated.hh"
#include "RegionFunctionColor.hh"
#include "RegionHalfSpace.hh"
#include "RegionLabeledSet.hh"
#include "RegionLineSegment.hh"
#include "RegionLogical.hh"
#include "RegionPlane.hh"
#include "RegionPoint.hh"
#include "RegionPolygon.hh"

#include "RegionFactory.hh"

namespace Amanzi {
namespace AmanziGeometry {

// Create region from XML specification
Teuchos::RCP<Region>
createRegion(const std::string reg_name,
             int reg_id,
             Teuchos::ParameterList& reg_spec,
             const Comm_type& comm)
{
  std::stringstream sstream;

  // There should be only one item below the region name which
  // indicates the shape of the region. Unfortunately, there is
  // nothing to prevent there from being multiple sublists.

  // There should be a "shape" or "region type" parameter instead of
  // this.  This should be fixed but would require a ton of change to
  // existing input files and the translator.  --etc

  if (reg_spec.numParams() != 1) {
    Errors::Message msg;
    msg << "Region spec \"" << reg_name
        << "\" should have exactly one shape sublist.";
    Exceptions::amanzi_throw(msg);
  }
  
  Teuchos::ParameterList::ConstIterator k = reg_spec.begin();
  std::string shape = reg_spec.name(k);
  std::string lifecycle_str = reg_spec.get<std::string>("Lifecycle", "Permanent");

  LifeCycleType lifecycle;
  if (lifecycle_str == "Permanent") {
    lifecycle = PERMANENT;
  } else if (lifecycle_str == "Temporary") {
    lifecycle = TEMPORARY;
  } else {
    Errors::Message msg;
    msg << "Region spec \"" << reg_name
        << "\" lifecycle must be either \"Permanent\" or \"Temporary\"";
    Exceptions::amanzi_throw(msg);
  }

  Teuchos::ParameterList& plist = reg_spec.sublist(shape);

  Teuchos::RCP<Region> region;
  if (shape == "region: box") {
    auto p0_vec = plist.get<Teuchos::Array<double> >("low coordinate");
    auto p1_vec = plist.get<Teuchos::Array<double> >("high coordinate");

    int dim = p0_vec.size();
    Point p0, p1;
      
    p0.set(dim, p0_vec.data());
    p1.set(dim, p1_vec.data());

    region = Teuchos::rcp(new RegionBox(reg_name, reg_id,
                                        p0, p1, lifecycle));

  } else if (shape == "region: plane") {
    auto p_vec = plist.get<Teuchos::Array<double> >("point");
    auto n_vec = plist.get<Teuchos::Array<double> >("normal");

    int dim = p_vec.size();
    Point p, n;
    p.set(dim, &(p_vec[0]));
    n.set(dim, &(n_vec[0]));

    region = Teuchos::rcp(new RegionPlane(reg_name, reg_id, p, n,
                                          lifecycle));

  } else if (shape == "region: halfspace") {
    auto p_vec = plist.get<Teuchos::Array<double> >("point");
    auto n_vec = plist.get<Teuchos::Array<double> >("normal");

    int dim = p_vec.size();
    Point p, n;
    p.set(dim, &(p_vec[0]));
    n.set(dim, &(n_vec[0]));

    region = Teuchos::rcp(new RegionHalfSpace(reg_name, reg_id, p, n,
                                              lifecycle));

  } else if (shape == "region: polygon") {
    int num_points = plist.get<int>("number of points");
    auto pvec = plist.get<Teuchos::Array<double> >("points");

    if (pvec.size()%num_points != 0) {
      Errors::Message mesg;
      mesg << "Incorrect number of values specified for polygon region \""
           << reg_name << "\"";        
      Exceptions::amanzi_throw(mesg);
    }

    int dim = pvec.size()/num_points;
    std::vector<Point> points;
    Point pnt(dim);
    for (int i = 0; i < num_points; ++i) {
      if (dim == 2)
        pnt.set(pvec[i*dim],pvec[i*dim+1]);
      else if (dim == 3)
        pnt.set(pvec[i*dim],pvec[i*dim+1],pvec[i*dim+2]);
      points.push_back(pnt);
    }

    region = Teuchos::rcp(new RegionPolygon(reg_name, reg_id,
                                            points, lifecycle));

  } else if (shape == "region: cylinder") {
    auto p_vec = plist.get<Teuchos::Array<double> >("point");
    auto a_vec = plist.get<Teuchos::Array<double> >("axis");
    double rad = plist.get<double>("radius");

    int dim = p_vec.size();
    Point p, a;
    p.set(dim, &(p_vec[0]));
    a.set(dim, &(a_vec[0]));

    region = Teuchos::rcp(new RegionCylinder(reg_name, reg_id, a, p, rad,
                                             lifecycle));

  } else if (shape == "region: labeled set") {
    std::string file = plist.get<std::string>("file");
    std::string format = plist.get<std::string>("format");
    std::string name = plist.get<std::string>("label");
    std::string entity_str = plist.get<std::string>("entity");

    if (entity_str == "Cell" || entity_str == "cell" || entity_str == "CELL")
      entity_str = "CELL";
    else if (entity_str == "Face" || entity_str == "face" || entity_str == "FACE")
      entity_str = "FACE";
    else if (entity_str == "Node" || entity_str == "node" || entity_str == "NODE")
      entity_str = "NODE";

    region = Teuchos::rcp(new RegionLabeledSet(reg_name, reg_id, entity_str,
                                               file, format, name, lifecycle));

  } else if (shape == "region: color function") {
    std::string file = plist.get<std::string>("file");
    int value = plist.get<int>("value");

    region = Teuchos::rcp(new RegionFunctionColor(reg_name, reg_id, file, 
                                                  value, comm, lifecycle));

  } else if (shape == "region: point") {
    auto p_vec = plist.get<Teuchos::Array<double> >("coordinate");
        
    int dim = p_vec.size();
    Point pnt;
    pnt.set(dim, &(p_vec[0]));

    region = Teuchos::rcp(new RegionPoint(reg_name, reg_id, pnt, lifecycle));

  } else if (shape == "region: logical") {
    std::string opstr = plist.get<std::string>("operation");
    auto region_names = plist.get< Teuchos::Array<std::string> >("regions");

    region = Teuchos::rcp(new RegionLogical(reg_name, reg_id, opstr,
                                            region_names.toVector(),
                                            lifecycle));

  } else if (shape == "region: enumerated set") {
    std::string entity_str = plist.get<std::string>("entity");
    std::vector<int> gids;

    if (plist.isParameter("file")) {
      std::string filename = plist.get<std::string>("file");
      HDF5Reader reader(filename);
      reader.ReadData(reg_name, gids);
    } else {
      auto entity_list = plist.get< Teuchos::Array<int> >("entity gids");
      gids = entity_list.toVector();
    }

    region = Teuchos::rcp(new RegionEnumerated(reg_name,
            reg_id, entity_str, gids, lifecycle));

  } else if (shape == "region: enumerated set from file") {
    Teuchos::ParameterList enum_params = reg_spec.sublist(shape);

    std::string filename = enum_params.get<std::string>("read from file");
    Teuchos::ParameterList enum_list_from_file = *Teuchos::getParametersFromXmlFile(filename);
    std::string enum_reg_name = enum_params.get<std::string>("region name");
    Teuchos::ParameterList enum_sub_list = enum_list_from_file.sublist(enum_reg_name);
    
    std::string entity_str = enum_sub_list.get<std::string>("entity");
    Teuchos::Array<int> entity_list =
      enum_sub_list.get< Teuchos::Array<int> >("entity gids");
    region = Teuchos::rcp(new RegionEnumerated(reg_name, reg_id, entity_str,
                                               entity_list.toVector(),
                                               lifecycle));
  } else if (shape == "region: all") {
    region = Teuchos::rcp(new RegionAll(reg_name, reg_id, lifecycle));

  } else if (shape == "region: boundary") {
    region = Teuchos::rcp(new RegionBoundary(reg_name, reg_id, lifecycle));

  } else if (shape == "region: box volume fractions") {
    auto p0_vec = plist.get<Teuchos::Array<double> >("corner coordinate");
    auto p1_vec = plist.get<Teuchos::Array<double> >("opposite corner coordinate");

    int dim = p0_vec.size();
    Teuchos::Array<double> normals_vec;

    if (plist.isParameter("normals")) {
      normals_vec = plist.get<Teuchos::Array<double> >("normals");
    } else {
      for (int i = 0; i < dim; ++i) {
        for (int j = 0; j < dim; ++j) {
          double tmp = (i == j) ? 1.0 : 0.0;
          normals_vec.push_back(tmp);
        }
      }
    }
    AMANZI_ASSERT(normals_vec.size() == dim*dim);

    Point p0, p1, p2;
    std::vector<Point> normals;
      
    p0.set(dim, p0_vec.data());
    p1.set(dim, p1_vec.data());
    for (int i = 0; i < dim; ++i) {
      p2.set(dim, &(normals_vec[i*dim]));
      normals.push_back(p2);
    }

    region = Teuchos::rcp(new RegionBoxVolumeFractions(
        reg_name, reg_id, p0, p1, normals, lifecycle));

  } else if (shape == "region: line segment") {
    auto p0_vec = plist.get<Teuchos::Array<double> >("end coordinate");
    auto p1_vec = plist.get<Teuchos::Array<double> >("opposite end coordinate");

    int dim = p0_vec.size();

    Point p0, p1;
    p0.set(dim, p0_vec.data());
    p1.set(dim, p1_vec.data());

    region = Teuchos::rcp(new RegionLineSegment(reg_name, reg_id, p0, p1, lifecycle));

  } else {
    Errors::Message mesg;
    mesg << "Cannot process region with shape \""
         << shape << "\"";
    Exceptions::amanzi_throw(mesg);
  }

  // tolerance for geometric operations
  double tolerance = TOL;
  if (plist.isSublist("expert parameters")) {
    tolerance = plist.sublist("expert parameters").get<double>("tolerance");
  }
  region->set_tolerance(tolerance);

  return region;
}


} // namespace AmanziGeometry
} // namespace Amanzi
