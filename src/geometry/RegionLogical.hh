/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Rao Garimella
*/

//! RegionLogical: A region defined by a logical operation on one or two other
//! regions

/*!

The list *region: logical* defines logical operations on regions allow for
more advanced region definitions. At this time the logical region allows for
logical operations on a list of regions.  *union* and *intersection* are
self-evident. In the case of *subtraction*, subtraction is performed from the
first region in the list.  The *complement* is a special case in that it is
the only case that operates on single region, and returns the complement to it
within the domain ENTIRE_DOMAIN.  Currently, multi-region booleans are not
supported in the same expression.

* `"operation`" ``[string]`` defines operation on the list of regions.
  Available options are *union*, *intersect*, *subtract*, *complement*

* `"regions`" ``[Array(string)]`` specifies the list of involved regions.

Example:

.. code-block:: xml

  <ParameterList name="LOWER_LAYERs">
    <ParameterList name="region: logical">
      <Parameter name="operation" type="string" value="union"/>
      <Parameter name="regions" type="Array(string)" value="{Middle1, Middle2,
Bottom}"/>
    </ParameterList>
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
  RegionLogical(const std::string& name, const int id,
                const std::string& operation_str,
                const std::vector<std::string>& component_regions,
                const LifeCycleType lifecycle = PERMANENT);

  // Label in the file
  BoolOpType operation() const { return operation_; }

  /// Is the the specified point inside this region
  bool inside(const Point& p) const;

  std::vector<std::string> component_regions() const
  {
    return component_regions_;
  }

 protected:
  BoolOpType operation_; // what logical operation should be performed
  std::vector<std::string> component_regions_; // names of regions in operation
};

} // namespace AmanziGeometry
} // namespace Amanzi

#endif
