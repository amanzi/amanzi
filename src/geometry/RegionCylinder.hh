/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
//! RegionHalfSpace: A planar (infinite) region in space, defined by a point and a normal.

/*
  Copyright 2010-2013 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Rao Garimella
           Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*!
List *region: cylinder* defines an infinite cylinder determined by 
a symmetry axis, point on this axis and radius.

* `"axis`" ``[Array(double)]`` symmetry axis

* `"point`" ``[Array(double)]`` point on a symmetry axis

* `"radius`" ``[double]`` cylinder radius

Example:

.. code-block:: xml

   <ParameterList name="TOP_SECTION"> <!-- parent list -->
     <ParameterList name="region: cylinder">
       <Parameter name="normal" type="Array(double)" value="{1, 1, 0}"/>
       <Parameter name="axis" type="Array(double)" value="{0, 0, 1}"/>
       <ParameterList name="expert parameters">
         <Parameter name="tolerance" type="double" value="1.0e-05"/>
       </ParameterList>
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
