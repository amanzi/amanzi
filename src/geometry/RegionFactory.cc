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

#include "Region.hh"
#include "RegionBox.hh"
#include "RegionPlane.hh"
#include "RegionLabeledSet.hh"
#include "RegionColorFunction.hh"
#include "RegionPoint.hh"
#include "RegionLogical.hh"
#include "RegionPolygon.hh"
#include "RegionEnumerated.hh"
#include "RegionAll.hh"
#include "RegionBoundary.hh"
#include "RegionBoxVolumeFractions.hh"

#include "RegionFactory.hh"

namespace Amanzi {
namespace AmanziGeometry {

// Create region from XML specification
Teuchos::RCP<Region>
createRegion(const std::string reg_name,
             Set_ID reg_id,
             Teuchos::ParameterList& reg_spec,
             const Epetra_MpiComm *comm)
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

  std::string lifecycle_str =
    reg_spec.get<std::string>("Lifecycle", "Permanent");
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

  Teuchos::RCP<Region> region;
  if (shape == "Region: Box") {
    Teuchos::ParameterList& box_params = reg_spec.sublist(shape);

    Teuchos::Array<double> p0_vec =
      box_params.get<Teuchos::Array<double> >("low coordinate");
        
    Teuchos::Array<double> p1_vec =
      box_params.get<Teuchos::Array<double> >("high coordinate");

    int dim = p0_vec.size();
    Point p0, p1;
      
    p0.set(dim,&(p0_vec[0]));
    p1.set(dim,&(p1_vec[0]));

    region = Teuchos::rcp(new RegionBox(reg_name, reg_id,
                                        p0, p1, lifecycle));

  } else if (shape == "Region: Plane") {
    Teuchos::ParameterList plane_params = reg_spec.sublist(shape);

    Teuchos::Array<double> p_vec =
      plane_params.get< Teuchos::Array<double> >("point");
        
    Teuchos::Array<double> n_vec =
      plane_params.get< Teuchos::Array<double> >("normal");

    int dim = p_vec.size();
    Point p, n;
    p.set(dim,&(p_vec[0]));
    n.set(dim,&(n_vec[0]));

    region = Teuchos::rcp(new RegionPlane(reg_name, reg_id, p, n,
                                          lifecycle));

  } else if (shape == "Region: Polygon") {
    Teuchos::ParameterList poly_params = reg_spec.sublist(shape);

    int num_points = poly_params.get<int>("number of points");
    Teuchos::Array<double> pvec =
      poly_params.get<Teuchos::Array<double> >("points");

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

  } else if (shape == "Region: Labeled Set") {
    Teuchos::ParameterList labeledset_params = reg_spec.sublist(shape);

    std::string file = labeledset_params.get<std::string>("file");
    std::string format = labeledset_params.get<std::string>("format");
    std::string name = labeledset_params.get<std::string>("label");
    std::string entity_str = labeledset_params.get<std::string>("entity");

    if (entity_str == "Cell" || entity_str == "cell" || entity_str == "CELL")
      entity_str = "CELL";
    else if (entity_str == "Face" || entity_str == "face" || entity_str == "FACE")
      entity_str = "FACE";
    else if (entity_str == "Node" || entity_str == "node" || entity_str == "NODE")
      entity_str = "NODE";

    region = Teuchos::rcp(new RegionLabeledSet(reg_name, reg_id, entity_str,
                                               file, format, name, lifecycle));

  } else if (shape == "Region: Color Function") {
    Teuchos::ParameterList colorfunc_params = reg_spec.sublist(shape);

    std::string file = colorfunc_params.get<std::string>("file");
    int value = colorfunc_params.get<int>("value");

    region = Teuchos::rcp(new RegionColorFunction(reg_name, reg_id, file, 
                                                  value, comm, lifecycle));

  } else if (shape == "Region: Point") {
    Teuchos::ParameterList point_params = reg_spec.sublist(shape);

    Teuchos::Array<double> p_vec =
      point_params.get<Teuchos::Array<double> >("coordinate");
        
    int dim = p_vec.size();
    Point pnt;
    pnt.set(dim,&(p_vec[0]));

    region = Teuchos::rcp(new RegionPoint(reg_name, reg_id, pnt, lifecycle));

  } else if (shape == "Region: Logical") {
    Teuchos::ParameterList logical_params = reg_spec.sublist(shape);
    std::string opstr = logical_params.get<std::string>("operation");
    Teuchos::Array<std::string> region_names = 
      logical_params.get< Teuchos::Array<std::string> >("regions");

    region = Teuchos::rcp(new RegionLogical(reg_name, reg_id, opstr,
                                            region_names.toVector(),
                                            lifecycle));

  } else if (shape == "Region: Enumerated Set") {
    Teuchos::ParameterList enum_params = reg_spec.sublist(shape);
    std::string entity_str = enum_params.get<std::string>("entity");

    if (entity_str == "Cell" || entity_str == "cell" || entity_str == "CELL")
      entity_str = "CELL";
    else if (entity_str == "Face" || entity_str == "face" || entity_str == "FACE")
      entity_str = "FACE";
    else if (entity_str == "Node" || entity_str == "node" || entity_str == "NODE")
      entity_str = "NODE";
      
    Teuchos::Array<int> entity_list = 
      enum_params.get< Teuchos::Array<int> >("entity gids");

    region = Teuchos::rcp(new RegionEnumerated(reg_name, reg_id, entity_str,
                                               entity_list.toVector(),
                                               lifecycle));

  } else if (shape == "Region: All") {
    region = Teuchos::rcp(new RegionAll(reg_name, reg_id, lifecycle));

  } else if (shape == "Region: Boundary") {
    region = Teuchos::rcp(new RegionBoundary(reg_name, reg_id, lifecycle));

  } else if (shape == "Region: Box Volume Fractions") {
    Teuchos::ParameterList& box_params = reg_spec.sublist(shape);

    Teuchos::Array<double> p0_vec =
        box_params.get<Teuchos::Array<double> >("corner coordinate");
        
    Teuchos::Array<double> p1_vec =
        box_params.get<Teuchos::Array<double> >("opposite corner coordinate");

    Teuchos::Array<double> normals_vec =
        box_params.get<Teuchos::Array<double> >("normals");

    int dim = p0_vec.size();
    Point p0, p1, p2;
    std::vector<Point> normals;
      
    p0.set(dim,&(p0_vec[0]));
    p1.set(dim,&(p1_vec[0]));
    for (int i = 0; i < dim; ++i) {
      p2.set(dim,&(normals_vec[i*dim]));
      normals.push_back(p2);
    }

    region = Teuchos::rcp(new RegionBoxVolumeFractions(
        reg_name, reg_id, p0, p1, normals, lifecycle));

  } else {
    Errors::Message mesg;
    mesg << "Cannot process region with shape \""
         << shape << "\"";
    Exceptions::amanzi_throw(mesg);
  }

  return region;
}


} // namespace AmanziGeometry
} // namespace Amanzi
