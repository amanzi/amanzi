/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Nonlinear finite volume method.
*/

#include "nlfv.hh"
#include "WhetStoneDefs.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* A harmonic averaging point (HAP) is a unique point on a plane 
* (line in 2D) * seprating two materials where (a) continuity 
* conditions are satisfied for continuous piecewise linear pressure
* functions and (b) pressure value is a convex combination of two 
* neighboring cell-based pressures p_c1 and p_c2:
*
*   p = w p_c1 + (1-w) p_c2. 
*
* Input: face f, two cells sharing this face, and two co-normal 
*        vectors Tni = Ti * normal where fixed normal is used.
* Output: HAP p and weight w.
****************************************************************** */
void
NLFV::HarmonicAveragingPoint(int f,
                             int c1,
                             int c2,
                             const AmanziGeometry::Point& Tn1,
                             const AmanziGeometry::Point& Tn2,
                             AmanziGeometry::Point& p,
                             double& weight)
{
  int dir;

  const AmanziGeometry::Point& fm = mesh_->face_centroid(f);
  const AmanziGeometry::Point& normal = mesh_->face_normal(f, false, c1, &dir);

  const AmanziGeometry::Point& cm1 = mesh_->cell_centroid(c1);
  const AmanziGeometry::Point& cm2 = mesh_->cell_centroid(c2);

  double area = mesh_->face_area(f);
  double d1 = fabs(normal * (fm - cm1)) / area;
  double d2 = fabs(normal * (fm - cm2)) / area;

  double a2 = area * area;
  double t1 = fabs(normal * Tn1) / a2;
  double t2 = fabs(normal * Tn2) / a2;

  double det = t1 * d2 + t2 * d1;
  weight = t1 * d2 / det;

  double factor = dir * d1 * d2 / det / area;
  p = weight * cm1 + (1.0 - weight) * cm2 + factor * (Tn1 - Tn2);
}


/* ******************************************************************
* Decomposion: conormal = w1 * tau[i1] + w2 * tau[i2] + w3 * tau[i3],
* where the weights ws = (w1, w2, w3) are non-negative.
****************************************************************** */
int
NLFV::PositiveDecomposition(int id1,
                            const std::vector<AmanziGeometry::Point>& tau,
                            const AmanziGeometry::Point& conormal,
                            double* ws,
                            int* ids)
{
  int ierr(1);
  int d = conormal.dim();
  int ntau = tau.size();

  // default is the TPFA stencil
  double c1 = norm(tau[id1]);
  double c2 = norm(conormal);

  ids[0] = id1;
  ws[0] = c2 / c1;
  for (int k = 1; k < d; k++) {
    ids[k] = (id1 + k) % ntau;
    ws[k] = 0.0;
  }

  // Check that the first direction is sufficient.
  double cs = (conormal * tau[id1]) / (c1 * c2);
  if (fabs(cs - 1.0) < WHETSTONE_TOLERANCE_DECOMPOSITION) return 0;

  // Find the other directions
  Tensor T(d, 2);
  double det(0.0);

  if (d == 2) {
    for (int i = 0; i < ntau; i++) {
      if (i == id1) continue;

      T.SetColumn(0, tau[id1]);
      T.SetColumn(1, tau[i]);

      // We skip almost colinear pairs.
      c2 = norm(tau[i]);
      double tmp = fabs(T.Det()) / (c1 * c2);
      if (tmp < WHETSTONE_TOLERANCE_DECOMPOSITION) continue;

      T.Inverse();
      AmanziGeometry::Point p = T * conormal;

      // We look for the strongest pair of vectors and try to
      // avoid degenerate cases to improve robustness.
      if (p[0] >= 0.0 && p[1] >= -WHETSTONE_TOLERANCE_DECOMPOSITION) {
        if (tmp > det) {
          det = tmp;
          ws[0] = p[0];
          ws[1] = fabs(p[1]);
          ids[1] = i;
          ierr = 0;
        }
      }
    }
  } else if (d == 3) {
    for (int i = 0; i < ntau; i++) {
      if (i == id1) continue;

      c2 = norm(tau[i]);
      for (int j = i + 1; j < ntau; j++) {
        if (j == id1) continue;

        T.SetColumn(0, tau[id1]);
        T.SetColumn(1, tau[i]);
        T.SetColumn(2, tau[j]);

        // We skip almost colinear pairs.
        double c3 = norm(tau[j]);
        double tmp = fabs(T.Det()) / (c1 * c2 * c3);
        if (tmp < WHETSTONE_TOLERANCE_DECOMPOSITION) continue;

        T.Inverse();
        AmanziGeometry::Point p = T * conormal;

        // We look for the strongest triple of vectors and try to
        // avoid degenerate cases to improve robustness.
        if (p[0] >= 0.0 && p[1] >= -WHETSTONE_TOLERANCE_DECOMPOSITION &&
            p[2] >= -WHETSTONE_TOLERANCE_DECOMPOSITION) {
          if (tmp > det) {
            det = tmp;
            ws[0] = p[0];
            ws[1] = fabs(p[1]);
            ws[2] = fabs(p[2]);
            ids[1] = i;
            ids[2] = j;
            ierr = 0;
          }
        }
      }
    }
  }

  return ierr;
}

} // namespace WhetStone
} // namespace Amanzi
