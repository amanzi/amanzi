/*
  WhetStone, version 2.0
  Release name: naka-to.

  Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Mimetic finite difference method.
*/


#include "nlfv.hh"
#include "WhetStoneDefs.hh"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* A harmonic averaging point is a unique point on a plane (line in 2D)
* seprating two materials where (a) continuity conditions are 
* satisfied for continuous piecewise linear pressure functions and 
* (b) pressure value is a convex combination of two neighboring 
* cell-based pressures p_0 and p_1:
*
*   p = w p_0 + (1-w) p_1. 
*
* NOTE: weigth is defined as 1.0 for a boundary face.
****************************************************************** */
void NLFV::HarmonicAveragingPoint(int face, std::vector<Tensor>& T,
                                  AmanziGeometry::Point& p, double& weight)
{
  int d = mesh_->space_dimension();

  Entity_ID_List cells;
  mesh_->face_get_cells(face, (ParallelTypeCast)WhetStone::USED, &cells);
  int ncells = cells.size();

  if (ncells == 1) {
    p = mesh_->face_centroid(face);
    weight = 1.0;
  } else {
    const AmanziGeometry::Point& fm = mesh_->face_centroid(face);
    const AmanziGeometry::Point& normal = mesh_->face_normal(face);

    const AmanziGeometry::Point& cm1 = mesh_->cell_centroid(cells[0]);
    const AmanziGeometry::Point& cm2 = mesh_->cell_centroid(cells[1]);

    AmanziGeometry::Point Tn1(d), Tn2(d);
    Tn1 = T[cells[0]] * normal;
    Tn2 = T[cells[1]] * normal;

    double d1 = fabs(normal * (fm - cm1));
    double d2 = fabs(normal * (fm - cm2));
    double t1 = fabs(normal * Tn1);
    double t2 = fabs(normal * Tn2);

    double det = t1 * d2 + t2 * d1;
    weight = t1 * d2 / det;

    AmanziGeometry::Point v1(d), v2(d);
    double area = mesh_->face_area(face);
    v1 = area * Tn1 - t1 * normal / area;
    v2 = area * Tn2 - t2 * normal / area;
    p = weight * cm1 + (1 - weight) * cm2 + (d1 * d2 / det) * (v2 - v1);
  }
}


/* ******************************************************************
* Decomposion: conormal = w1 * tau[i1] + w2 * tau[i2] + w3 * tau[i3],
* where the weights ws = (w1, w2, w3) are non-negative.
****************************************************************** */
int NLFV::PositiveDecomposition(
    int id1, const std::vector<AmanziGeometry::Point>& tau,
    const AmanziGeometry::Point& conormal, double* ws, int* ids)
{
  int d = mesh_->space_dimension();
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
  double det = 0.0;

  if (d == 2) {
    for (int i = 0; i < ntau; i++) {
      if (i == id1) continue;

      T.SetRow(0, tau[id1]);
      T.SetRow(1, tau[i]);

      // We skip almost colinear pairs.
      c2 = norm(tau[i]);
      double tmp = fabs(T.Det()) / (c1 * c2);
      if (tmp < WHETSTONE_TOLERANCE_DECOMPOSITION) continue;

      T.Inverse();
      AmanziGeometry::Point p = T * conormal;

      // We look for the strongest pair of vectors and try to 
      // avoid degenerate cases to improve robustness.
      if (p[0] >= 0.0 && 
          p[1] >= -WHETSTONE_TOLERANCE_DECOMPOSITION) {
        if (tmp > det) {
          det = tmp;
          ws[0] = p[0];
          ws[1] = fabs(p[1]);
          ids[1] = i; 
        }
      }
    }
  } else if (d == 3) {
    for (int i = 0; i < ntau; i++) {
      if (i == id1) continue;

      c2 = norm(tau[i]);
      for (int j = i + 1; j < ntau; j++) {
        if (j == id1) continue;

        T.SetRow(0, tau[id1]);
        T.SetRow(1, tau[i]);
        T.SetRow(2, tau[j]);

        // We skip almost colinear pairs.
        double c3 = norm(tau[j]);
        double tmp = fabs(T.Det()) / (c1 * c2 * c3);
        if (tmp < WHETSTONE_TOLERANCE_DECOMPOSITION) continue;

        T.Inverse();
        AmanziGeometry::Point p = T * conormal;

        // We look for the strongest triple of vectors and try to 
        // avoid degenerate cases to improve robustness.
        if (p[0] >= 0.0 && 
            p[1] >= -WHETSTONE_TOLERANCE_DECOMPOSITION &&
            p[2] >= -WHETSTONE_TOLERANCE_DECOMPOSITION) {
          if (tmp > det) {
            det = tmp;
            ws[0] = p[0];
            ws[1] = fabs(p[1]);
            ws[2] = fabs(p[2]);
            ids[1] = i; 
            ids[2] = j; 
          }
        }
      }
    }
  }

  return 0;
}

}  // namespace WhetStone
}  // namespace Amanzi

