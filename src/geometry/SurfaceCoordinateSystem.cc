/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Useful tools for local surface coordinate system.
*/

#include <cmath>

#include "SurfaceCoordinateSystem.hh"

namespace Amanzi {
namespace AmanziGeometry {

/* ******************************************************************
* Initialization of local coordinate system
****************************************************************** */
void
SurfaceCoordinateSystem::Init()
{
  int d = normal_.dim();
  tau_ = std::make_shared<std::vector<AmanziGeometry::Point>>(d - 1);
  auto& tau = *tau_;

  normal_unit_ /= AmanziGeometry::norm(normal_unit_);

  if (d == 2) {
    tau[0] = AmanziGeometry::Point(-normal_[1], normal_[0]);
  } else {
    // here we need orthogonal projector to implement hierarchical construction
    if (fabs(normal_unit_[0]) > 0.1)
      tau[0] = AmanziGeometry::Point(normal_unit_[1], -normal_unit_[0], 0.0);
    else
      tau[0] = AmanziGeometry::Point(0.0, -normal_unit_[2], normal_unit_[1]);

    tau[0] /= norm(tau[0]);
    tau[1] = normal_unit_ ^ tau[0];
  }
}


/* ******************************************************************
* Initialization of local coordinate system
****************************************************************** */
AmanziGeometry::Point
SurfaceCoordinateSystem::Project(const AmanziGeometry::Point& x, bool flag) const
{
  int d = tau_->size();
  AmanziGeometry::Point xloc(d);

  for (int i = 0; i < d; ++i) xloc[i] = (flag) ? ((x - origin_) * (*tau_)[i]) : (x * (*tau_)[i]);

  return xloc;
}

} // namespace AmanziGeometry
} // namespace Amanzi
