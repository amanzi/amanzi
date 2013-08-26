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
                                  AmanziGeometry::Point& hap, double& hap_weight)
{
  int d = mesh_->space_dimension();

  AmanziMesh::Entity_ID_List cells;
  mesh_->face_get_cells(face, AmanziMesh::USED, &cells);
  int ncells = cells.size();

  if (ncells == 1) {
    hap = mesh_->face_centroid(face);
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
    hap_weight = t1 * d2 / det;

    AmanziGeometry::Point v1(d), v2(d);
    double area = mesh_->face_area(face);
    v1 = area * Tn1 - t1 * normal / area;
    v2 = area * Tn2 - t2 * normal / area;
    hap = hap_weight * cm1 + (1 - hap_weight) * cm2 + (d1 * d2 / det) * (v2 - v1);
  }
}

}  // namespace WhetStone
}  // namespace Amanzi

