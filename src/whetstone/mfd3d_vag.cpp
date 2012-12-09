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

#include <cmath>
#include <vector>

#include "mfd3d.hpp"
#include "tensor.hpp"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* A harmonic point is a unique point on a plane seprating two 
* materials where (a) continuity conditions are satisfied for 
* continuous piecewise linear pressure functions and (b) pressure 
* value is a convex combination of two neighboring cell-based 
* prossures.
****************************************************************** */
void MFD3D::CalculateHarmonicPoints(int face,
                                    std::vector<Tensor>& T,
                                    AmanziGeometry::Point& harmonic_point,
                                    double& harmonic_point_weight)
{
  int d = mesh_->space_dimension();

  AmanziMesh::Entity_ID_List cells;
  mesh_->face_get_cells(face, AmanziMesh::USED, &cells);
  int ncells = cells.size();

  if (ncells == 1) {
    harmonic_point = mesh_->face_centroid(face);
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
    harmonic_point_weight = t1 * d2 / det;

    AmanziGeometry::Point v1(d), v2(d);
    double area = mesh_->face_area(face);
    v1 = area * Tn1 - t1 * normal / area;
    v2 = area * Tn2 - t2 * normal / area;
    harmonic_point = harmonic_point_weight * cm1
                   + (1 - harmonic_point_weight) * cm2
                   + (d1 * d2 / det) * (v2 - v1);
  }
}


/* ******************************************************************
* Dispesion flux is based on values at first d harmonic points. 
* The flux is scaled by the area of corresponding subface (facet).
* The case of Dirichlet-Neumann corner has to be implemented 
* separately using array bc_face_id().(lipnikov@lanl.gov).
****************************************************************** */
int MFD3D::DispersionCornerFluxes(int node, int cell, Tensor& dispersion,
                                  std::vector<AmanziGeometry::Point>& corner_points,
                                  double cell_value,
                                  std::vector<double>& corner_values,
                                  std::vector<double>& corner_fluxes)
{
  int d = mesh_->space_dimension();

  AmanziMesh::Entity_ID_List faces, nodes;
  mesh_->node_get_cell_faces(node, cell, AmanziMesh::USED, &faces);
  int nfaces = faces.size();
  const AmanziGeometry::Point& cm = mesh_->cell_centroid(cell);

  // gradient calculation grad(C) = X^{-1} * d_values
  Tensor X(d, 2);
  AmanziGeometry::Point dp[d];

  for (int i = 0; i < d; i++) {
    (dp[i]).init(d);
    dp[i] = corner_points[i] - cm;
    for (int j = 0; j < d; j++) X(j, i) = (dp[i])[j];
  }
  X.inverse();

  AmanziGeometry::Point gradient(d), dvalues(d);
  for (int i = 0; i < d; i++) dvalues[i] = corner_values[i] - cell_value;
  gradient = X * dvalues;
  gradient = dispersion * gradient;

  corner_fluxes.clear();
  for (int i = 0; i < nfaces; i++) {  // calculate corner fluxes
    int f = faces[i];
    mesh_->face_get_nodes(f, &nodes);
    int nnodes = nodes.size();

    const AmanziGeometry::Point& normal = mesh_->face_normal(f);
    corner_fluxes.push_back((gradient * normal) / nnodes);
  }

  return 0;
}

}  // namespace WhetStone
}  // namespace Amanzi
 
