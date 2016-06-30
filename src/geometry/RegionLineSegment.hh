/*

 A region is defined by line segment, or in other words by two point in space.
 It includes all cells that intersect with the line segment

 Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatsky (dasvyat@lanl.gov)

*/

#ifndef AMANZI_LINE_SEGMENT_REGION_HH
#define AMANZI_LINE_SEGMENT_REGION_HH

#include <vector>

#include "Point.hh"
#include "Region.hh"

namespace Amanzi {
namespace AmanziGeometry {

  class RegionLineSegment: public Region {
  public:
// Default constructor uses two corner points (order not important)
    RegionLineSegment(const std::string& name,
                      const Set_ID id,
                      const Point& p0, 
                      const Point& p1,
                      const LifeCycleType lifecycle=PERMANENT);

// Is the specified point inside this region?
  bool inside(const Point& p) const;

  // Calculate intersection polytope object with the line segment.
  // Return value is 1 if the line segment intersects a polytope,
  // a 0 if not.
  //
  // Polyhedron with counter clockwise ordered faces (wrt normals)
  double intersect(const std::vector<Point>& polytope, 
                   const std::vector<std::vector<int> >& faces) const;

 protected:
    const Point p0_, p1_; // two end points of the line.
    std::vector<Point> line_points_;
    std::vector<double> line_frac_;
    std::vector<int> reg_entities_;

  };

double PlaneLineIntersection(const std::vector<Point>& plane,
                             const std::vector<Point>& line);
  
double det_aux(const std::vector<double>& first_row,
               const std::vector<double>& submatr);

} // namespace AmanziGeometry
} // namespace Amanzi

#endif
