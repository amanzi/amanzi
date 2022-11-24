/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/*
  A region defined by a set of mesh entities in a mesh file

  Strictly speaking, we should tie this region class to a particular
  mesh or mesh file but that cause a circular dependency of meshes
  on regions and of labeled set regions on meshes. We will rely on the 
  fact that when a mesh is created specifying a geometric model, it
  will create mesh entity sets based on the labeled sets in that 
  geometric model. 

  If we need to change this behavior, then we can make a forward
  declaration of AmanziMesh::Mesh, make the Mesh class a friend, add
  a mesh variable to this class and have a protected method to set
  the mesh

  The region will consist of all mesh elements for which the indicator
  function is a particular value at their centroids

  Copyright 2010-2013 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Rao Garimella
*/

#include "dbc.hh"
#include "errors.hh"

#include "RegionLabeledSet.hh"

namespace Amanzi {
namespace AmanziGeometry {

// -------------------------------------------------------------
// RegionLabeledSet:: constructor
// -------------------------------------------------------------
RegionLabeledSet::RegionLabeledSet(const std::string& name,
                                   const int id,
                                   const std::string& entity_str,
                                   const std::string& file,
                                   const std::string& format,
                                   const std::string& label,
                                   const LifeCycleType lifecycle)
  : Region(name, id, false, RegionType::LABELEDSET, 0, 0, lifecycle),
    entity_str_(entity_str),
    file_(file),
    format_(format),
    label_(label)
{}


// -------------------------------------------------------------
// RegionLabeledSet::inside
// -------------------------------------------------------------
bool
RegionLabeledSet::inside(const Point& p) const
{
  Errors::Message mesg("In/out check not implemented for labeled sets");
  Exceptions::amanzi_throw(mesg);
  return false;
}

} // namespace AmanziGeometry
} // namespace Amanzi
