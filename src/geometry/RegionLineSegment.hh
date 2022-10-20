/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
//! RegionLineSegment: A line segment, defined by two points in space.

/*
 Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatsky (dasvyat@lanl.gov)
*/

/*!

List *region: line segment* desribes a region defined by a line
segment. This region is a set of cells which intersect with a line
segment.  The line segment is allowed to intersect with one or more cells. Zero length
line segments are allowed. The line segment is defined by its ends
points.

.. _region-line-segment-spec:
.. admonition:: region-line-segment-spec

    * `"end coordinate`" ``[Array(double)]`` Location of one end of a line
      segment.
    * `"opposite end coordinate`" ``[Array(double)]`` Location of the opposite
      end of a line segment.

Example:

.. code-block:: xml

   <ParameterList name="WELL"> <!-- parent list -->
      <ParameterList name="region: line segment">
        <Parameter name="end coordinate" type="Array(double)" value="{497542.44, 5393755.77, 0.0}"/>
        <Parameter name="opposite end coordinate" type="Array(double)" value="{497542.44, 5393755.77, 100.0}"/>
      </ParameterList>
    </ParameterList>     

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
                    const int id,
                    const Point& p0, 
                    const Point& p1,
                    const LifeCycleType lifecycle=LifeCycleType::PERMANENT);

  // Is the specified point inside this region?
  bool inside(const Point& p) const;

  // Calculate intersection polytope object with the line segment.
  // Return value is 1 if the line segment intersects a polytope,
  // a 0 if not.
  //
  // Polyhedron with counter clockwise ordered faces (wrt normals)
  double intersect(const std::vector<Point>& polytope, 
                   const std::vector<std::vector<int> >& faces) const;

  void ComputeInterLinePoints(const std::vector<Point>& polytope,
                              const std::vector<std::vector<int> >& faces,
                              Point& res_point) const;

 protected:
  const Point p0_, p1_; // two end points of the line.
  //std::vector<Point> line_points_;
  //std::vector<double> line_frac_;
  //bool complete_;
};

double PlaneLineIntersection(const std::vector<Point>& plane,
                             const std::vector<Point>& line);
  
double det_aux(const std::vector<double>& first_row,
               const std::vector<double>& submatr);

} // namespace AmanziGeometry
} // namespace Amanzi

#endif
