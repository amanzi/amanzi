/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
//! RegionBox: a rectangular region in space, defined by two points

/*
  Copyright 2010-2013 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: William Perkins
           Rao Garimella
           Ethan Coon (ecoon@lanl.gov)

*/

/*!

List *region: box* defines a region bounded by coordinate-aligned
planes. Boxes are allowed to be of zero thickness in only one
direction in which case they are equivalent to planes.

* `"low coordinate`" ``[Array(double)]`` Location of the boundary point with the lowest coordinates.

* `"high coordinate`" ``[Array(double)]`` Location of the boundary points with the highest coordinates.

Example:

.. code-block:: xml

   <ParameterList name="WELL">  <!-- parent list -->
     <ParameterList name="region: box">
       <Parameter name="low coordinate" type="Array(double)" value="{-5.0,-5.0, -5.0}"/>
       <Parameter name="high coordinate" type="Array(double)" value="{5.0, 5.0,  5.0}"/>
     </ParameterList>
   </ParameterList>
  
*/


#ifndef AMANZI_BOX_REGION_HH_
#define AMANZI_BOX_REGION_HH_

#include "Point.hh"

#include "Region.hh"

namespace Amanzi {
namespace AmanziGeometry {

class RegionBox : public Region {
 public:

  // Default constructor uses two corner points (order not important).
  RegionBox(const std::string& name,
            const Set_ID id,
            const Point& p0, 
            const Point& p1,
            const LifeCycleType lifecycle=PERMANENT);

  // Get the first point defining the region
  const Point& point0() const { return p0_; }

  // Get the second point defining the region
  const Point& point1() const { return p1_; }

  // Is the the specified point inside this region
  bool inside(const Point& p) const;

  // corners
  void corners(Point *lo_corner, Point *hi_corner) const;

  // Is the box degenerate - zero length in one or more directions and
  // if so in how many directions?
  bool is_degenerate(int *ndeg) const;

protected:
  
  const Point p0_; // one corner of the region
  const Point p1_; // the other corner of the region

  // Is the specified value between the two values (inclusive, order not important)
  //
  // static -- this could be a function, but is only used here and
  // therefore scoping is convenient.
  static bool between_(const double& x, const double& x0, const double& x1);

};


} // namespace AmanziGeometry
} // namespace Amanzi


#endif
