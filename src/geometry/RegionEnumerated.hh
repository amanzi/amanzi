/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon (coonet@ornl.gov)
*/

//! RegionEnumerated: A region enumerated as a list of IDs.

/*!

List *region: enumerated set* defines a set of mesh entities via the list
of input global ids. Note that global ids are not defined correctly when
parallle mesh is created on a fly.

* `"entity`" ``[string]`` Type of the mesh object.  Valid are *cell*, *face*,
*edge*, *node*

* `"entity gids`" ``[Array(int)]`` List of the global IDs of the entities.


Example:

.. code-block:: xml

   <ParameterList name="WELL"> <!-- parent list -->
     <ParameterList name="region: enumerated set">
       <Parameter name="entity" type="string" value="face"/>
       <Parameter name="entity gids" type="Array(int)" value="{1, 12, 23, 34}"/>
     </ParameterList>
   </ParameterList>

*/

#ifndef AMANZI_ENUMERATED_REGION_HH_
#define AMANZI_ENUMERATED_REGION_HH_

#include <vector>

#include "errors.hh"
#include "GeometryDefs.hh"

#include "Region.hh"

namespace Amanzi {
namespace AmanziGeometry {

class RegionEnumerated : public Region {
 public:
  RegionEnumerated(const std::string& name, const int id,
                   const std::string& entity_str,
                   const std::vector<Entity_ID>& ents,
                   const LifeCycleType lifecycle = PERMANENT);

  const std::vector<Entity_ID>& entities() const { return entities_; }
  const std::string& entity_str() const { return entity_str_; }

  // Is the the specified point inside this region
  bool inside(const Point& p) const;

 protected:
  std::string entity_str_; // what kind of entities make up this set
  const std::vector<Entity_ID> entities_; // list of those included
};


} // namespace AmanziGeometry
} // namespace Amanzi


#endif
