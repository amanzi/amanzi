/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/*
  A point in space

  Copyright 2010-2013 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Rao Garimella
*/

#ifndef AMANZI_REGION_POINT_HH_
#define AMANZI_REGION_POINT_HH_

#include "Region.hh"

namespace Amanzi {
namespace AmanziGeometry {

class RegionPoint : public Region {
 public:
  RegionPoint(const std::string& name,
              const Set_ID id,
              const Point& p,
              const LifeCycleType lifecycle=PERMANENT);

  // Get the point defining the region
  const Point& point(void) const { return p_; }

  // Is the specified point inside this region - in this case it
  // means coincident with the region point
  bool inside(const Point& p) const;

 protected:
  const Point p_;
};

} // namespace AmanziGeometry
} // namespace Amanzi


#endif
