/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! An enumerated set is simply a list of discrete mesh global IDs.
/*!

An enumerated set and a labeled set are functionally equivalent, but where a
labeled set stores its labels in the mesh file (typically as an attribute on
the mesh), an enumerated set stores its labeles separately, either directly in
the input file (for short lists) or in a separate file (for long lists).

Note that global ids are not defined correctly when a parallel mesh is created on
a fly.  Instead, prepartition the mesh and use the partitioned global IDs.

`"region type`" = `"enumerated set`"

.. _region-enumerated-set-spec:
.. admonition:: region-enumerated-set-spec

   * `"entity`" ``[string]`` Type of the mesh object.  One of: `"cell`",
     `"face`", `"edge`", `"node`"
   * `"entity gids`" ``[Array(int)]`` List of the global IDs of the entities.


Example:

.. code-block:: xml

   <ParameterList name="WELL"> <!-- parent list -->
     <Parameter name="region type" type="string" value="enumerated set"/>
     <Parameter name="entity" type="string" value="face"/>
     <Parameter name="entity gids" type="Array(int)" value="{1, 12, 23, 34}"/>
   </ParameterList>


For an enumerated set from file, the file should be an XML file that contains a
sublist whose name is the region name.  That sublist should be a
`"region-enumerated-set-spec`".

`"region type`" = `"enumerated set from file`"
   
.. _region-enumerated-set-from-file-spec:
.. admonition:: region-enumerated-set-from-file-spec

   * `"filename`" ``[string]`` Name of the XML file to read.

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
  RegionEnumerated(const std::string& name,
                   const int id,
                   const std::string& entity_str,
                   const std::vector<Entity_ID>& ents,
                   const LifeCycleType lifecycle = LifeCycleType::PERMANENT);

  const std::vector<Entity_ID>& entities() const { return entities_; }
  const std::string& entity_str() const { return entity_str_; }

  // Is the the specified point inside this region
  bool inside(const Point& p) const;

 protected:
  std::string entity_str_;                // what kind of entities make up this set
  const std::vector<Entity_ID> entities_; // list of those included
};


} // namespace AmanziGeometry
} // namespace Amanzi


#endif
