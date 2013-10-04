/*
  This is the mimetic discretization component of the Amanzi code. 

  Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Release name: aka-to.
  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/


#include "WhetStoneDefs.hh"
#include "nlfv.hh"


namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* A harmonic point is a unique point on a plane seprating two 
* materials where (a) continuity conditions are satisfied for 
* continuous piecewise linear pressure functions and (b) pressure 
* value is a convex combination of two neighboring cell-based 
* prossures.
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
    weight = 0.5;
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
* where the weights w1,w2,w3 >= 0.
* We make three assumptions:
* 1. The directions tau_i are normalized, ||tau_i||=1.
* 2. The first direction is fized.
* 3. conormal can be rewritten. 
****************************************************************** */
int NLFV::PositiveDecomposition(
    int id1, const std::vector<AmanziGeometry::Point>& tau,
    AmanziGeometry::Point& conormal, double* ws, int* ids)
{
  int d = mesh_->space_dimension();
  int ntau = tau.size();

  double cnorm = norm(conormal);
  conormal /= cnorm;

  // elliminate first direction
  ids[0] = id1;

  double w1 = conormal * tau[id1];  
  ws[0] = w1 * cnorm;
  if (w1 < 0.0) return -1;  // Perhaps this cell is non-convex.
  for (int k = 0; k < d; k++) conormal[k] -= w1 * tau[id1][k];
  
  if (norm(conormal) < WHETSTONE_TOLERANCE_DECOMPOSITION) {
    for (int k = 1; k < d; k++) {
      ids[k] = (id1 + k) % ntau;
      ws[k] = 0.0;
    }
    return 0;
  }

  // find the second direction
  int mtau = ntau - 1;
  int jds[mtau];

  int id2 = -1;
  double w2 = 0.0;
  for (int i = 0; i < mtau; i++) {
    int k = (id1 + i) % ntau;
    jds[i] = k;
    double tmp = conormal * tau[k];
    if (tmp > w2) {
      w2 = tmp;
      id2 = k; 
    }
  }
  ids[1] = id2;
  ws[1] = w2;
  if (id2 < 0) return -1;  // Perhaps this cell is non-convex.
  for (int k = 0; k < d; k++) conormal[k] -= w2 * tau[id2][k];

  // find the third direction
  if (d == 3) {
    if (norm(conormal) < WHETSTONE_TOLERANCE_DECOMPOSITION) {
      ids[2] = jds[(id2 + 1) % mtau];
      ws[2] = 0.0;
      return 0;
    }
 
    ntau--;
    mtau--;

    int id3 = -1;
    double w3 = 0.0;
    for (int i = 0; i < mtau; i++) {
      int k = (id2 + i) % ntau;
      double tmp = conormal * tau[k];
      if (tmp > w3) {
        w3 = tmp;
        id3 = k; 
      }
    }
    ids[2] = id3;
    ws[2] = w3;
    if (id3 < 0) return -1;  // Perhaps this cell is non-convex.
  }

  return 0;
}

}  // namespace WhetStone
}  // namespace Amanzi

