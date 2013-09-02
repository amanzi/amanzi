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
****************************************************************** */
void NLFV::PositiveDecomposition(const AmanziGeometry::Point& conormal, 
                                 const std::vector<AmanziGeometry::Point>& tau,
                                 double* w1, double* w2, int* i1, int* i2)
{
  int n = tau.size();

  // calculate non-negative projections
  std::vector<int> prj_id;
  std::vector<double> prj_w;
  
  for (int i = 0; i < n; i++) {
    double a = conormal * tau[i];
    if (a >= 0.0) {
      prj_id.push_back(i);
      prj_w.push_back(i);
    }
  }

  // find extrema for these projections
  *i1 = prj_id[0];
  *i2 = *i1;

  int m = prj_id.size();
  for (int i = 1; i < m; i++) {
  }

  *w1 = prj_w[*i1];
  *w2 = prj_w[*i2];
}


}  // namespace WhetStone
}  // namespace Amanzi

