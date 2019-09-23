/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Useful tools for various local coordinate systems.
*/

#ifndef AMANZI_COORDINATE_SYSTEMS_HH_
#define AMANZI_COORDINATE_SYSTEMS_HH_

#include <cmath>

#include "Point.hh"

namespace Amanzi {
namespace WhetStone {

// -- comparison operators
inline void FaceCoordinateSystem(const AmanziGeometry::Point& normal,
                                 std::vector<AmanziGeometry::Point>& tau) 
{
  int d = normal.dim();
  tau.resize(d - 1);    
  if (d == 2) {
    tau[0] = AmanziGeometry::Point(-normal[1], normal[0]);
  } else {
    double area = norm(normal);
    if (fabs(normal[0]) > 0.1)
      tau[0] = AmanziGeometry::Point(normal[1], -normal[0], 0.0);
    else 
      tau[0] = AmanziGeometry::Point(0.0, -normal[2], normal[1]);

    tau[0] *= sqrt(area) / norm(tau[0]);
    tau[1] = normal ^ tau[0];
    tau[1] /= area;
  }
}

}  // namespace WhetStone
}  // namespace Amanzi

#endif
