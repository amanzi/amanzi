/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Useful tools for a local surface coordinate system: local coordinate
  system and transformation of coodinates to and from this system.
*/

#include <cmath>

#include "SurfaceCoordinateSystem.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* Initialization of local coordinate system
****************************************************************** */
void SurfaceCoordinateSystem::Init()
{
  int d = normal_.dim();
  tau_ = std::make_shared<std::vector<AmanziGeometry::Point> >(d - 1);
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
* Project a vector on manifold.
****************************************************************** */
AmanziGeometry::Point SurfaceCoordinateSystem::Project(
   const AmanziGeometry::Point& x, bool flag) const
{
  int d = tau_->size();
  AmanziGeometry::Point xloc(d);    

  for (int i = 0; i < d; ++i) 
    xloc[i] = (flag) ? ((x - origin_) * (*tau_)[i]) : (x * (*tau_)[i]);

  return xloc;
}


/* ******************************************************************
* Inverse transformation. We could use the SVD but direct algorithm
* could be faster for small values of d.
****************************************************************** */
DenseMatrix SurfaceCoordinateSystem::InverseTransformation() const
{
  const auto& B = *tau_;
 
  int d = B.size();
  int dnew = B[0].dim();

  DenseMatrix Binv(dnew, d);
  Binv.PutScalar(0.0);

  if (d == 2) {
    // find monor with the largest determinant
    int i0(0), i1(1);
    double det01, det02, det12, tmp01, tmp02, tmp12;
    
    det01 = B[0][0] * B[1][1] - B[0][1] * B[1][0];
    det02 = B[0][0] * B[1][2] - B[0][2] * B[1][0];
    det12 = B[0][1] * B[1][2] - B[0][2] * B[1][1];

    tmp01 = std::fabs(det01); 
    tmp02 = std::fabs(det02);
    tmp12 = std::fabs(det12);

    if (tmp02 >= std::max(tmp01, tmp12)) {
      i1 = 2;
      det01 = det02;
    } else if (tmp12 >= std::max(tmp01, tmp02)) {
      i0 = 1;
      i1 = 2;
      det01 = det12;
    }

    Binv(i0, 0) = B[1][i1] / det01;
    Binv(i1, 0) =-B[1][i0] / det01;

    Binv(i0, 1) =-B[0][i1] / det01;
    Binv(i1, 1) = B[0][i0] / det01;
  }

  return Binv;
}

}  // namespace WhetStone
}  // namespace Amanzi

