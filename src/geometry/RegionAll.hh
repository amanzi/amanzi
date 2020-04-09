/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon (coonet@ornl.gov)
*/

//! A region consisting of all entities on a mesh.

#ifndef AMANZI_REGION_ALL_HH_
#define AMANZI_REGION_ALL_HH_

#include <vector>

#include "errors.hh"
#include "GeometryDefs.hh"

#include "Region.hh"

namespace Amanzi {
namespace AmanziGeometry {

class RegionAll : public Region {
 public:
  RegionAll(const std::string& name, const int id,
            const LifeCycleType lifecycle = PERMANENT);

  // Is the the specified point inside this region
  bool inside(const Point& p) const;
};

} // namespace AmanziGeometry
} // namespace Amanzi


#endif
