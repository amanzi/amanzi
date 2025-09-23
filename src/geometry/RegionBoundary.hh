/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

//! A region consisting of the entire boundary of any mesh.
/*!

No parameters are required.

`"region type`" = `"boundary`"

.. _region-boundary-spec:
.. admonition:: region-boundary-spec

   * `"empty`" ``[bool]`` **True** This is simply here to avoid issues with empty lists.

Example:

.. code-block:: xml

   <ParameterList name="DOMAIN_BOUNDARY"> <!-- parent list names the region -->
     <Parameter name="region type" type="string" value="boundary"/>
   </ParameterList>

*/

#ifndef AMANZI_REGION_BOUNDARY_HH_
#define AMANZI_REGION_BOUNDARY_HH_

#include "errors.hh"
#include "GeometryDefs.hh"
#include "Region.hh"

namespace Amanzi {
namespace AmanziGeometry {

class RegionBoundary : public Region {
 public:
  RegionBoundary(const std::string& name,
                 const int id,
                 const LifeCycleType lifecycle = LifeCycleType::PERMANENT);

  // Is the the specified point inside this region
  bool inside(const Point& p) const;
};

} // namespace AmanziGeometry
} // namespace Amanzi

#endif
