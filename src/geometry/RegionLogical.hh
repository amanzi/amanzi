/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Rao Garimella
*/

//! A region defined by a logical operation on one or more other regions.
/*!

The logical region defines logical operations on regions to allow for more
advanced region definitions. *union* and *intersection* are self-evident. In
the case of *subtraction*, all elements in the first region are included, then
all elements in every other region after the first are removed.  The
*complement* returns all entries that are not in the given region.

In all cases, logical operations are done discretely on the resolved
regions. For instance, the union of two box regions is the discrete union of
each box, resolved on a mesh, not the the geometric union of the two boxes.


`"region type`" = `"logical`"

.. _region-logical-spec:
.. admonition:: region-logical-spec

   * `"operation`" ``[string]`` defines operation on the list of regions.
     One of: `"union`", `"intersect`", `"subtract`", `"complement`"
   * `"regions`" ``[Array(string)]`` specifies the list of involved regions.

Example:

.. code-block:: xml

   <ParameterList name="LOWER_LAYERs">
     <Parameter name="region type" type="string" value="logical"/>
     <Parameter name="operation" type="string" value="union"/>
     <Parameter name="regions" type="Array(string)" value="{Middle1, Middle2, Bottom}"/>
   </ParameterList>

*/

#ifndef AMANZI_REGION_LOGICAL_HH_
#define AMANZI_REGION_LOGICAL_HH_

#include "Region.hh"

namespace Amanzi {
namespace AmanziGeometry {

class RegionLogical : public Region {
 public:
  // constructor
  RegionLogical(const std::string& name,
                const int id,
                const std::string& operation_str,
                const std::vector<std::string>& component_regions,
                const LifeCycleType lifecycle = LifeCycleType::PERMANENT);

  // Label in the file
  BoolOpType get_operation() const { return operation_; }

  /// Is the the specified point inside this region
  bool inside(const Point& p) const;

  const std::vector<std::string>& get_component_regions() const { return component_regions_; }

 protected:
  BoolOpType operation_;                       // what logical operation should be performed
  std::vector<std::string> component_regions_; // names of regions in operation
};

} // namespace AmanziGeometry
} // namespace Amanzi

#endif
