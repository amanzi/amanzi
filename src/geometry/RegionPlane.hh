/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/*
  A planar (infinite) region in space, defined by a point and a plane

  Copyright 2010-2013 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Rao Garimella
           Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_REGION_PLANE_HH_
#define AMANZI_REGION_PLANE_HH_

#include "Region.hh"

namespace Amanzi {
namespace AmanziGeometry {

class RegionPlane : public Region {

 public:
  // Default constructor uses point and normal
  RegionPlane(const std::string& name,
              const Set_ID id,
              const Point& p, 
              const Point& normal,
              const LifeCycleType lifecycle=PERMANENT);
  
  // Get the point defining the plane
  const Point& point() const { return p_; }

  // Get the normal point defining the plane
  const Point& normal() const { return n_; }

  // Is the specified point inside this region - in this case it
  // means on the plane
  bool inside(const Point& p) const;

protected:
  const Point p_;              /* point on the plane */
  const Point n_;              /* normal to the plane */

};

} // namespace AmanziGeometry
} // namespace Amanzi


#endif
