/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Rao Garimella
*/

//! A single point in space.
/*!

A point region defines a point in space.  Note that this region will match all
cells that contain this point in the closure of that cell.

Note that the dimension of the coordinate must match that of the mesh on which
the region is resolved.

`"region type`" = `"point`"

.. _region-point-spec:
.. admonition:: region-point-spec

   * `"coordinate`" ``[Array(double)]`` Location of point in space.

Example:

.. code-block:: xml

   <ParameterList name="DOWN_WIND150">
     <Parameter name="region type" type="string" value="point"/>
     <Parameter name="coordinate" type="Array(double)" value="{-150.0, 0.0, 0.0}"/>
   </ParameterList>

*/

#ifndef AMANZI_REGION_POINT_HH_
#define AMANZI_REGION_POINT_HH_

#include "Region.hh"

namespace Amanzi {
namespace AmanziGeometry {

class RegionPoint : public Region {
 public:
  RegionPoint(const std::string& name,
              const int id,
              const Point& p,
              const LifeCycleType lifecycle = LifeCycleType::PERMANENT);

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
