/*
This is the mimetic discretization component of the Amanzi code. 

Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided Reconstruction.cppin the top-level COPYRIGHT file.

Release name: aka-to.
Author: Konstantin Lipnikov (lipnikov@lanl.gov)
Usage: 
*/


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
* Decompose: conormal = w1 * tau[i1] + w2 * tau[i2],  w1,w2 >= 0.
* Requires some modification.
****************************************************************** */
void NLFV::MaximumDecomposition(const AmanziGeometry::Point& conormal, 
                                const std::vector<AmanziGeometry::Point>& tau,
                                double* w1, double* w2, int* i1, int* i2)
{
  // calculate all projections
  int ntau = tau.size();
  double cs[ntau];
  
  for (int i = 0; i < ntau; i++) {
    cs[i] = conormal * tau[i];
  }

  // find the largest value in these projections
  int k = 0;
  for (int i = 1; i < ntau; i++) {
    if (cs[i] > cs[k]) k = i;
  }
  *i1 = k;
  *w1 = cs[k];

  // find the closest to zero in these projections
  k = 0;
  for (int i = 1; i < ntau; i++) {
    if (fabs(cs[i]) < fabs(cs[k])) k = i;
  }
  *i2 = k;
  *w2 = cs[k];
}

}  // namespace WhetStone
}  // namespace Amanzi

