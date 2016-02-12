/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/*
  Abstract class for a Region, which is a geometric or discrete
  subdomain, along with some basic setters/getters.

  A geometric region is just some arbitrary subset of space, that can
  be specified in a myriad of ways.  At a minimum, there is a need to
  be able to determine if a point is inside that space.  A disrete
  region is an enumerated (via input spec or labeled sets inside the
  mesh file) list of entities, and is specific to a mesh.

  The region class does not use a constructor based on the XML parameter
  list because it has to create derived region classes based on the shape 
  parameter of the region specification.

  Copyright 2010-2013 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: William Perkins
           Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_REGION_HH_
#define AMANZI_REGION_HH_

#include <string>

#include "VerboseObject.hh"
#include "GeometryDefs.hh"

namespace Amanzi {
namespace AmanziGeometry {

class Point;

class Region {
 public:

  // Destructor
  virtual ~Region() {}

  // Dimension of the subdomain
  unsigned int manifold_dimension() const {
    return manifold_dimension_;
  }
  void set_manifold_dimension(unsigned int dimension) {
    manifold_dimension_ = dimension;
  }

  // Dimension of points in the subdomain
  unsigned int space_dimension() const {
    return space_dimension_;
  }
  void set_space_dimension(unsigned int dimension) {
    space_dimension_ = dimension;
  }
  
  // Name of the region -- no setter (set by constructor)
  std::string name() const {
    return name_;
  }

  // Integer identifier of the region
  Set_ID id() const {
    return id_;
  }
  void set_id(Set_ID id) {
    id_ = id;
  }

  // Geometric/enumerated
  bool is_geometric() const {
    return geometric_;
  }
  
  // Type of the region
  RegionType type() const {
    return type_;
  }

  // Get the Lifecycle of this region - Do mesh entity sets derived from
  // it have to be kept around or are they temporary and can be destroyed
  // as soon as they are used?
  LifeCycleType lifecycle() const 
  {
    return lifecycle_;
  }

  // Is the specified point inside the closure of the Region
  virtual bool inside(const Point& p) const = 0;

 protected:

  // Constructor -- protected as it should never be called directly
  Region(const std::string& name,
         Set_ID id,
         bool geometric,
         RegionType type,
         unsigned int dim,
         unsigned int geom_dim,
         LifeCycleType lifecycle=PERMANENT)
    : name_(name),
      id_(id),
      geometric_(geometric),
      type_(type),
      manifold_dimension_(dim),
      space_dimension_(geom_dim),
      lifecycle_(lifecycle) {}
    

 protected:

  // Lifecycle (Temporary or Permanent)
  LifeCycleType lifecycle_;
  
  // Topological dimension of region (0, 1, 2, 3)
  unsigned int manifold_dimension_;
  unsigned int space_dimension_;

  // Name of identifier
  std::string name_;

  // Integer identifier of region
  Set_ID id_;

  // Region type
  RegionType type_;

  // Geometric or enumerated?
  bool geometric_;

 private:
  Region(const Region& other); // prevent copy constructor
  Region& operator=(const Region& other); // prevent operator=
  
  
};



typedef std::vector<Region> RegionVector;

} // namespace AmanziGeometry
} // namespace Amanzi

#endif

