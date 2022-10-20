/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
//! RegionHalfSpace: A planar (infinite) region in space, defined by a point and a normal.

/*
  Copyright 2010-2013 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Rao Garimella
           Ethan Coon (ecoon@lanl.gov)
*/

/*!
List *region: halfspace* defines a halfspace determined by a plane and 
outward normal to the plane.

* `"normal`" ``[Array(double)]`` Normal to the plane.

* `"point`" ``[Array(double)]`` Point in space.

Example:

.. code-block:: xml

   <ParameterList name="TOP_SECTION"> <!-- parent list -->
     <ParameterList name="region: halfspace">
       <Parameter name="point" type="Array(double)" value="{2, 3, 5}"/>
       <Parameter name="normal" type="Array(double)" value="{1, 1, 0}"/>
       <ParameterList name="expert parameters">
         <Parameter name="tolerance" type="double" value="1.0e-05"/>
       </ParameterList>
     </ParameterList>
   </ParameterList>

*/


#ifndef AMANZI_REGION_HALFSPACE_HH_
#define AMANZI_REGION_HALFSPACE_HH_

#include "Region.hh"

namespace Amanzi {
namespace AmanziGeometry {

class RegionHalfSpace : public Region {

 public:
  // Default constructor uses point and normal
  RegionHalfSpace(const std::string& name,
                  const int id,
                  const Point& p, 
                  const Point& normal,
                  const LifeCycleType lifecycle=LifeCycleType::PERMANENT);
  
  // Get the point defining the plane
  const Point& point() const { return p_; }

  // Get the normal point defining the plane
  const Point& normal() const { return n_; }

  // Is the specified point inside this region - in this case it
  // means on the plane
  bool inside(const Point& p) const;

protected:
  const Point p_;  /* point on the plane */
  const Point n_;  /* normal to the plane */
};

} // namespace AmanziGeometry
} // namespace Amanzi


#endif
