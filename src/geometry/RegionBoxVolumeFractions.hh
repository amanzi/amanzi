/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Lipnikov Konstantin (lipnikov@lanl.gov)
      Rao Garimella (rao@lanl.gov)
*/

//! RegionBoxVolumeFractions: A rectangular region in space, defined by two corner points and normals to sides.
/*!

List *region: box volume fraction* defines a region bounded by a box *not*
aligned with coordinate axes.
Boxes are allowed to be of zero thickness in only one direction in which case
they are equivalent to rectangles on a plane or segments on a line.

.. _region-box-volume-fractions-spec:
.. admonition:: region-box-volume-fractions-spec

    * `"corner coordinate`" ``[Array(double)]`` Location of one box corner.
    * `"opposite corner coordinate`" ``[Array(double)]`` Location of the opposite box corner.
    * `"normals`" ``[Array(double)]`` Normals to sides in a linear array. Default is columns of
      the identity matrix. The normals may be scaled arbitrarily but must be orthogonal to
      one another and form the right coordinate frame.

Example:

.. code-block:: xml

   <ParameterList name="BASIN">  <!-- parent list -->
     <ParameterList name="region: box volume fractions">
       <Parameter name="corner coordinate" type="Array(double)" value="{-1.0,-1.0, 1.0}"/>
       <Parameter name="opposite corner coordinate" type="Array(double)" value="{1.0, 1.0, 1.0}"/>
       <Parameter name="normals" type="Array(double)" value="{1.0, 0.0, 0.0
                                                              0.0, 2.0, 0.0,
                                                              0.0, 0.0, 3.0}"/>
     </ParameterList>
   </ParameterList>

This example defines a degenerate box, a square on a surface *z=1*.

 */

#ifndef AMANZI_BOX_VOLUME_FRACTIONS_REGION_HH_
#define AMANZI_BOX_VOLUME_FRACTIONS_REGION_HH_

#include <vector>
#include <list>

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
                           const int id,
                           const Point& p0,
                           const Point& p1,
                           const std::vector<Point>& normals,
                           const LifeCycleType lifecycle = LifeCycleType::PERMANENT);

  // Is the specified point inside this region?
  bool inside(const Point& p) const;

  // Calculate intersection volume of polytope object with the box. The operation
  // is well defined when the polytope and box have same dimensionality.
  //
  // Polyhedron with counter clockwise ordered faces (wrt normals)
  double
  intersect(const std::vector<Point>& polytope, const std::vector<std::vector<int>>& faces) const;

 protected:
  const Point p0_, p1_; // two corners of the box
  const std::vector<Point> normals_;

 private:
  TensorSimple N_;
  double jacobian_; // change of area/volume during transformation
  int degeneracy_;  // degenerate box direction, only is allowed.
};


struct ClippedFace {
  ClippedFace()
  {
    new_edge = std::make_pair(-1, -1);
    edge_flag = 0;
  }

  std::list<int> nodes;
  std::pair<int, int> new_edge;
  int edge_flag;
};


// non-member functions
// -- intersection of counter clockwise oriented convex polygons
void
IntersectConvexPolygons(const std::vector<Point>& xy1,
                        const std::vector<Point>& xy2,
                        std::vector<Point>& xy3);

// -- intersection of a convex polyhedra, one is defined by a set of half-spaces
void
IntersectConvexPolyhedra(const std::vector<Point>& xyz1,
                         const std::vector<std::vector<int>>& faces1,
                         const std::vector<std::pair<Point, Point>>& xyz2,
                         std::vector<Point>& xyz3,
                         std::vector<std::vector<int>>& faces3);

} // namespace AmanziGeometry
} // namespace Amanzi

#endif
