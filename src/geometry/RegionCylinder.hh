/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Rao Garimella
           Konstantin Lipnikov (lipnikov@lanl.gov)
*/

//! A cylindrical region in 3 dimensions.
/*!

A cylinder region defines an infinite cylinder determined by
a symmetry axis, point on this axis, and radius.

`"region type`" = `"cylinder`"

.. _region-cylinder-spec:
.. admonition:: region-cylinder-spec

   * `"axis`" ``[Array(double)]`` symmetry axis
   * `"point`" ``[Array(double)]`` point on a symmetry axis
   * `"radius`" ``[double]`` cylinder radius

Example:

.. code-block:: xml

   <ParameterList name="TOP_SECTION"> <!-- parent list -->
     <Parameter name="region type" type="string" value="cylinder"/>
     <Parameter name="normal" type="Array(double)" value="{1, 1, 0}"/>
     <Parameter name="axis" type="Array(double)" value="{0, 0, 1}"/>
     <ParameterList name="expert parameters">
       <Parameter name="tolerance" type="double" value="1.0e-05"/>
     </ParameterList>
   </ParameterList>

*/


#ifndef AMANZI_REGION_CYLINDER_HH_
#define AMANZI_REGION_CYLINDER_HH_

#include "Region.hh"

namespace Amanzi {
namespace AmanziGeometry {

class RegionCylinder : public Region {
 public:
  // Default constructor uses point and normal
  RegionCylinder(const std::string& name,
                 const int id,
                 const Point& axis,
                 const Point& p,
                 double radius,
                 const LifeCycleType lifecycle = LifeCycleType::PERMANENT);

  // Is the specified point inside this region - in this case it
  // means on the plane
  bool inside(const Point& p) const;

 protected:
  const Point p_;    // point on the plane
  const Point axis_; // axis of symmetry
  double rad2_;
};

} // namespace AmanziGeometry
} // namespace Amanzi


#endif
