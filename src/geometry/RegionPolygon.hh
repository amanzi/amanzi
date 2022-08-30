/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
//! RegionPolygon: A closed polygonal segment of a plane.

/*
  Copyright 2010-2013 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Rao Garimella
*/

/*!

The list *region: polygon* defines a polygonal region on which mesh faces and
nodes can be queried. NOTE that one cannot ask for cells in a polygonal surface
region. In 2D, the polygonal region is a line and is specified by 2 points.
In 3D, the polygonal region is specified by an arbitrary number of points.
In both cases the point coordinates are given as a linear array. The polygon
can be non-convex.

This provides a set of faces with a normal for computing flux.

The polygonal surface region can be queried for a normal. In 2D, the normal is
defined as [Vy,-Vx] where [Vx,Vy] is the vector from point 1 to point 2.
In 3D, the normal of the polygon is defined by the order in which points 
are specified.

``[region-polygon-spec]``

* `"number of points`" ``[int]`` Number of polygon points.

* `"points`" ``[Array(double)]`` Point coordinates in a linear array. 

Example:

.. code-block:: xml

   <ParameterList name="XY_PENTAGON">
     <ParameterList name="region: polygon">
       <Parameter name="number of points" type="int" value="5"/>
       <Parameter name="points" type="Array(double)" value="{-0.5, -0.5, -0.5, 
                                                              0.5, -0.5, -0.5,
                                                              0.8, 0.0, 0.0,
                                                              0.5,  0.5, 0.5,
                                                             -0.5, 0.5, 0.5}"/>
       <ParameterList name="expert parameters">
         <Parameter name="tolerance" type="double" value="1.0e-3"/>
       </ParameterList>
     </ParameterList>
   </ParameterList>

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
                const int id, 
                const std::vector<Point>& polypoints, 
                const LifeCycleType lifecycle=LifeCycleType::PERMANENT);

  typedef std::vector<AmanziGeometry::Point>::const_iterator PointIterator;
  std::size_t size() const { return points_.size(); }
  PointIterator begin() const { return points_.begin(); }
  PointIterator end() const { return points_.end(); }

  const Point& normal() const { return normal_; }

  // Is the the specified point inside this region
  bool inside(const Point& p) const;

protected:
  
  std::vector<Point> points_;  /* Points of the polygon */
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
