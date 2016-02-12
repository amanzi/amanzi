/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/*
  A closed polygonal segment of a plane.

  Copyright 2010-2013 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Rao Garimella
*/

#ifndef AMANZI_REGION_POLYGON_HH_
#define AMANZI_REGION_POLYGON_HH_

#include "Region.hh"

namespace Amanzi {
namespace AmanziGeometry {

class RegionPolygon : public Region {
 public:
  // Default constructor uses two corner points (order not important).
  RegionPolygon(const std::string& name,
                const Set_ID id, 
                const std::vector<Point>& polypoints, 
                const LifeCycleType lifecycle=PERMANENT);
  
  unsigned int num_points() const { return points_.size(); }
  const std::vector<const Point>& points() const { return points_; }

  const Point& normal() const { return normal_; }

  // Is the the specified point inside this region
  bool inside(const Point& p) const;

protected:
  
  std::vector<const Point> points_;  /* Points of the polygon */
  Point normal_;                     /* Normal to the polygon */
  unsigned int elim_dir_;            /* Coord dir to eliminate while projecting
                                        polygon for in/out tests 
                                        0 - yz, eliminate x coord        
                                        1 - xz, eliminate y coord        
                                        2 - xy, eliminate z coord        */

private:
  void Init_();
};

} // namespace AmanziGeometry
} // namespace Amanzi


#endif
