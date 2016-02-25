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

#ifndef AMANZI_REGION_LABELED_SET_HH_
#define AMANZI_REGION_LABELED_SET_HH_

#include "Region.hh"

namespace Amanzi {
namespace AmanziGeometry {

class RegionLabeledSet : public Region {
 public:
  // constructor 
  RegionLabeledSet(const std::string& name, 
                   const Set_ID id, 
                   const std::string& entity_str,
                   const std::string& file,
                   const std::string& format,
                   const std::string& label,
                   const LifeCycleType lifecycle=PERMANENT);
  
  // Label in the file
  const std::string& label() const { return label_; }

  // Is the the specified point inside this region
  bool inside(const Point& p) const;

  const std::string& entity_str() const { return entity_str_; }

protected:  
  const std::string entity_str_; // what kind of entities make up this set
  const std::string file_; // which file are we supposed to read it from
  const std::string format_; // format of the file
  const std::string label_; // Label used to identify set in the file (may be different from name)
};

} // namespace AmanziGeometry
} // namespace Amanzi


#endif
