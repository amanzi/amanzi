/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
//! RegionPoint: a point in space.

/*
  Copyright 2010-2013 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Rao Garimella
*/


/*!
List *region: point* defines a point in space. 
This region consists of cells containing this point.

.. _region-point-spec:
.. admonition:: region-point-spec

    * `"coordinate`" ``[Array(double)]`` Location of point in space.

Example:

.. code-block:: xml

   <ParameterList name="DOWN_WIND150"> <!-- parent list defining the name -->
     <ParameterList name="region: point">
       <Parameter name="coordinate" type="Array(double)" value="{-150.0, 0.0, 0.0}"/>
     </ParameterList>
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
              const LifeCycleType lifecycle=PERMANENT);

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
