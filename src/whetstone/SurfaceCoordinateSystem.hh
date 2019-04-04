/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Useful tools for local surface coordinate system.
*/

#ifndef AMANZI_SURFACE_COORDINATE_SYSTEM_HH_
#define AMANZI_SURFACE_COORDINATE_SYSTEM_HH_

#include "Point.hh"

namespace Amanzi {
namespace WhetStone {

// -- comparison operators
class SurfaceCoordinateSystem {
 public:
  SurfaceCoordinateSystem(const AmanziGeometry::Point& normal)
    : normal_(normal) { Init(); }
  ~SurfaceCoordinateSystem() {};

  // calculate orthogonal vectors of the surface coordinate system
  void Init();

  // project vector on the surface
  AmanziGeometry::Point Project(const AmanziGeometry::Point& x) const;

  // const access
  const AmanziGeometry::Point& normal() const { return normal_; }
  const std::shared_ptr<std::vector<AmanziGeometry::Point> > tau() const { return tau_; }

 private:
  AmanziGeometry::Point normal_;
  std::shared_ptr<std::vector<AmanziGeometry::Point> > tau_;
};

}  // namespace WhetStone
}  // namespace Amanzi

#endif
