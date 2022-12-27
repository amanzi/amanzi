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

#ifndef AMANZI_SURFACE_COORDINATE_SYSTEM_HH_
#define AMANZI_SURFACE_COORDINATE_SYSTEM_HH_

#include <memory>

#include "Point.hh"

namespace Amanzi {
namespace WhetStone {

// -- comparison operators
class SurfaceCoordinateSystem {
 public:
  SurfaceCoordinateSystem(const AmanziGeometry::Point& origin, const AmanziGeometry::Point& normal)
    : origin_(origin), normal_(normal), normal_unit_(normal)
  {
    Init();
  }
  ~SurfaceCoordinateSystem(){};

  // calculate orthogonal vectors of the surface coordinate system
  void Init();

  // project vector on the surface.
  // -- flag = true, calculate coordinates of new points relative to origin
  // -- flag = false, project vectors
  AmanziGeometry::Point Project(const AmanziGeometry::Point& x, bool flag) const;

  // const access
  const AmanziGeometry::Point& get_origin() const { return origin_; }
  const AmanziGeometry::Point& normal() const { return normal_; }
  const AmanziGeometry::Point& normal_unit() const { return normal_unit_; }
  const std::shared_ptr<std::vector<AmanziGeometry::Point>> tau() const { return tau_; }

 private:
  AmanziGeometry::Point origin_, normal_, normal_unit_;
  std::shared_ptr<std::vector<AmanziGeometry::Point>> tau_;
};

} // namespace WhetStone
} // namespace Amanzi

#endif
