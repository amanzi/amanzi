/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
//! RegionBoundary:  A region consisting of all entities on the domain boundary

/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*!

List *region: boundary* defines a set of all boundary faces. 
Using this definition, faces located on the domain boundary are extracted.

.. _region-boundary-spec:
.. admonition:: region-boundary-spec

    * `"entity`" ``[string]`` Type of the mesh object.  Unclear whether this is
      used or can be other things than `"face`"?

Example:

.. code-block:: xml

   <ParameterList name="DOMAIN_BOUNDARY"> <!-- parent list names the region -->
     <ParameterList name="region: boundary">
       <Parameter name="entity" type="string" value="face"/>
     </ParameterList>
   </ParameterList>

*/

#ifndef AMANZI_REGION_BOUNDARY_HH_
#define AMANZI_REGION_BOUNDARY_HH_

#include <vector>

#include "errors.hh"
#include "GeometryDefs.hh"
#include "Region.hh"

namespace Amanzi {
namespace AmanziGeometry {

class RegionBoundary : public Region {
 public:
  RegionBoundary(const std::string& name,
                 const int id,
                 const LifeCycleType lifecycle=PERMANENT);

  // Is the the specified point inside this region
  bool inside(const Point& p) const;
  
 protected:
  const std::string entity_str_; // what kind of entities make up this set
  const std::vector<int> entities_; // list of those included
};

}  // namespace AmanziGeometry
}  // namespace Amanzi

#endif
