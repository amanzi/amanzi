/*
  A rectangular region in space, defined by two corner points and 
  normals to sides.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Lipnikov Konstantin (lipnikov@lanl.gov)
           Rao Garimella (rao@lanl.gov)
*/

#ifndef AMANZI_BOX_VOLUME_FRACTIONS_REGION_HH_
#define AMANZI_BOX_VOLUME_FRACTIONS_REGION_HH_

#include <vector>

#include "Point.hh"
#include "Region.hh"
#include "TensorSimple.hh"

namespace Amanzi {
namespace AmanziGeometry {

class RegionBoxVolumeFractions : public Region {
 public:
  // Default constructor uses two corner points (order not important)
  // and vector of normals (order is important).
  RegionBoxVolumeFractions(const std::string& name,
                           const Set_ID id,
                           const Point& p0, 
                           const Point& p1,
                           const std::vector<Point>& normals,
                           const LifeCycleType lifecycle=PERMANENT);

  // Is the specified point inside this region?
  bool inside(const Point& p) const;

  // Calculate intersection volume of polytope object with the box. The operation
  // is well defined when the polytope and box have same dimensionality.
  double intersect(const std::vector<Point>& polytope) const;

 protected:
  const Point p0_, p1_; // two corners of the box
  const std::vector<Point> normals_;

 private:
  TensorSimple N_;
  double jacobian_;  // change of area/volume during transformation
  int degeneracy_;  // degenerate box direction, only is allowed.
};


// non-member functions
// -- intersection of counter clockwise oriented convex polygons
void IntersectConvexPolygons(const std::vector<Point>& xy1,
                             const std::vector<Point>& xy2,
                             std::vector<Point>& xy3);

// -- intersection of a convex polyhedra, one is defined by a set of half-spaces
void IntersectConvexPolyhedra(const std::vector<Point>& xyz1,
                              const std::vector<std::vector<int> >& faces1,
                              const std::vector<std::pair<Point, Point> >& xyz2,
                              std::vector<Point>& xyz3,
                              std::vector<std::vector<int> >& faces3);

}  // namespace AmanziGeometry
}  // namespace Amanzi

#endif
