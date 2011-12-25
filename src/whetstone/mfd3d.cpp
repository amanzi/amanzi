/*
This is the mimetic discretization component of the Amanzi code. 
License: BSD
Release name: aka-to.
Author: Konstantin Lipnikov (lipnikov@lanl.gov)
Usage: 
*/

#include <cmath>

#include "Teuchos_SerialDenseMatrix.hpp"

#include "Mesh.hh"
#include "Point.hh"

#include "mfd3d.hpp"
#include "tensor.hpp"

namespace Amanzi {
namespace WhetStone {

/* ******************************************************************
* Darcy mass matrix: a wrapper for other low-level routines
****************************************************************** */
int MFD3D::darcy_mass(int cell,
                       Tensor& permeability,
                       Teuchos::SerialDenseMatrix<int, double>& M)
{
  int d = mesh_->space_dimension();

  AmanziMesh::Entity_ID_List faces;
  mesh_->cell_get_faces(cell, &faces);
  int nfaces = faces.size();
 
  Teuchos::SerialDenseMatrix<int, double> N(nfaces, d);
  Teuchos::SerialDenseMatrix<int, double> Mc(nfaces, nfaces);

  int ok = L2_consistency(cell, permeability, N, Mc);
  if (ok) return WHETSTONE_ELEMENTAL_MATRIX_WRONG;

  stability_scalar(cell, N, Mc, M);
  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* Darcy mass matrix: a wrapper for other low-level routines
****************************************************************** */
int MFD3D::darcy_mass_inverse(int cell,
                              Tensor& permeability,
                              Teuchos::SerialDenseMatrix<int, double>& W)
{
  int d = mesh_->space_dimension();

  AmanziMesh::Entity_ID_List faces;
  mesh_->cell_get_faces(cell, &faces);
  int nfaces = faces.size();
 
  Teuchos::SerialDenseMatrix<int, double> R(nfaces, d);
  Teuchos::SerialDenseMatrix<int, double> Wc(nfaces, nfaces);

  int ok = L2_consistency_inverse(cell, permeability, R, Wc);
  if (ok) return WHETSTONE_ELEMENTAL_MATRIX_WRONG;

  stability_scalar(cell, R, Wc, W);
  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* A harmonic point is a unique point on a plane seprating two 
* materials where (a) continuity conditions are satisfied for 
* continuous piecewise linear pressure functions and (b) pressure 
* value is a convex combination of two neighboring cell-based 
* prossures.
****************************************************************** */
void MFD3D::calculate_harmonic_points(int face, 
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
  }
  else {
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
int MFD3D::dispersion_corner_fluxes(int node,
                                    int cell,
                                    Tensor& dispersion,
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

  for (int i=0; i<d; i++) {
    (dp[i]).init(d);
    dp[i] = corner_points[i] - cm;
    for (int j=0; j<d; j++) X(j, i) = (dp[i])[j];  
  }
  X.inverse();

  AmanziGeometry::Point gradient(d), dvalues(d);
  for (int i=0; i<d; i++) dvalues[i] = corner_values[i] - cell_value;
  gradient = X * dvalues;
  gradient = dispersion * gradient;

  corner_fluxes.clear();
  for (int i=0; i<nfaces; i++) {  // calculate corner fluxes
    int f = faces[i];
    mesh_->face_get_nodes(f, &nodes);
    int nnodes = nodes.size();

    const AmanziGeometry::Point& normal = mesh_->face_normal(f);
    corner_fluxes.push_back((gradient * normal) / nnodes);
  }

  return 0;
}


/* ******************************************************************
* Consistency condition for inner product in space of Darcy fluxes. 
* Only upper triangular part of Wc is calculated.
****************************************************************** */
int MFD3D::L2_consistency(int cell,
                          Tensor& T,
                          Teuchos::SerialDenseMatrix<int, double>& N,
                          Teuchos::SerialDenseMatrix<int, double>& Mc)
{
  AmanziMesh::Entity_ID_List faces;
  mesh_->cell_get_faces(cell, &faces);
 
  int num_faces = faces.size();
  if (num_faces != N.numRows()) return num_faces;  // matrix was not reshaped

  int d = mesh_->space_dimension();
  double volume = mesh_->cell_volume(cell);

  AmanziGeometry::Point v1(d), v2(d);
  const AmanziGeometry::Point& cm = mesh_->cell_centroid(cell);
 
  Tensor Tinv(T);
  Tinv.inverse();
 
  for (int i=0; i<num_faces; i++) {
    int f = faces[i];
    const AmanziGeometry::Point& fm = mesh_->face_centroid(f);
    double a1 = mesh_->face_area(f);

    v1 = a1 * (fm - cm);
    v2 = Tinv * v1;
 
    for (int j=i; j<num_faces; j++) {
      f = faces[j];
      const AmanziGeometry::Point& fm = mesh_->face_centroid(f);
      double a2 = mesh_->face_area(f);

      v1 = a2 * (fm - cm);
      Mc(i, j) = (v1 * v2) / volume; 
    }
  }
 
  std::vector<int> dirs;
  mesh_->cell_get_face_dirs(cell, &dirs);

  for (int i=0; i<num_faces; i++) {
    int f = faces[i];
    const AmanziGeometry::Point& normal = mesh_->face_normal(f);
    double area = dirs[i] * mesh_->face_area(f);

    for (int k=0; k<d; k++) N(i, k) = normal[k] / area;
  }
  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* Consistency condition for inverse of mass matrix in space of Darcy 
* fluxes. Only the upper triangular part of Wc is calculated.
****************************************************************** */
int MFD3D::L2_consistency_inverse(int cell,
                                  Tensor& T,
                                  Teuchos::SerialDenseMatrix<int, double>& R,
                                  Teuchos::SerialDenseMatrix<int, double>& Wc)
{
  AmanziMesh::Entity_ID_List faces;
  std::vector<int> dirs;

  mesh_->cell_get_faces(cell, &faces);
  mesh_->cell_get_face_dirs(cell, &dirs);

  int num_faces = faces.size();
  if (num_faces != R.numRows()) return num_faces;  // matrix was not reshaped

  int d = mesh_->space_dimension();
  AmanziGeometry::Point v1(d);
  double volume = mesh_->cell_volume(cell);

  for (int i=0; i<num_faces; i++) {
    int f = faces[i];
    const AmanziGeometry::Point& normal = mesh_->face_normal(f);
    double a1 = dirs[i] * mesh_->face_area(f);

    v1 = T * normal;
 
    for (int j=i; j<num_faces; j++) {
      f = faces[j];
      const AmanziGeometry::Point& v2 = mesh_->face_normal(f);
      double a2 = dirs[j] * mesh_->face_area(f);

      Wc(i, j) = (v1 * v2) / (a1 * a2 * volume); 
    }
  }

  const AmanziGeometry::Point& cm = mesh_->cell_centroid(cell);

  for (int i=0; i<num_faces; i++) {
    int f = faces[i];
    const AmanziGeometry::Point& fm = mesh_->face_centroid(f);
    double area = mesh_->face_area(f);

    for (int k=0; k<d; k++) R(i, k) = area * (fm[k] - cm[k]);
  }
  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* Consistency condition for stifness matrix in heat conduction. 
* Only the upper triangular part of Ac is calculated.
****************************************************************** */
int MFD3D::H1_consistency(int cell,
                          Tensor& T,
                          Teuchos::SerialDenseMatrix<int, double>& N,
                          Teuchos::SerialDenseMatrix<int, double>& Ac)
{
  AmanziMesh::Entity_ID_List nodes, faces;

  mesh_->cell_get_nodes(cell, &nodes);
  int num_nodes = nodes.size();
  if (num_nodes != N.numRows()) return num_nodes;  // matrix was not reshaped

  mesh_->cell_get_faces(cell, &faces);
  int num_faces = faces.size();

  std::vector<int> dirs;
  mesh_->cell_get_face_dirs(cell, &dirs);

  int d = mesh_->space_dimension();
  double volume = mesh_->cell_volume(cell);
  AmanziGeometry::Point p(d), pnext(d), pprev(d), v1(d), v2(d), v3(d);

  /* to calculate matrix R, we use temporary matrix N */
  N = 0;

  for (int i=0; i<num_faces; i++) { 
    int f = faces[i];
    const AmanziGeometry::Point& normal = mesh_->face_normal(f);
    const AmanziGeometry::Point& fm = mesh_->face_centroid(f);
    double area = mesh_->face_area(f);
    
    AmanziMesh::Entity_ID_List face_nodes; 
    mesh_->face_get_nodes(f, &face_nodes);
    int num_face_nodes = face_nodes.size();

    for (int j=0; j<num_face_nodes; j++) {
      int jnext = (j + 1) % num_face_nodes;
      int jprev = (j + num_face_nodes - 1) % num_face_nodes;

      int v = face_nodes[j];
      int vnext = face_nodes[jnext];
      int vprev = face_nodes[jprev];

      mesh_->node_get_coordinates(v, &p);
      mesh_->node_get_coordinates(vnext, &pnext);
      mesh_->node_get_coordinates(vprev, &pprev);

      v1 = pprev - pnext;
      v2 = p - fm;
      v3 = v1^v2;
      double u = dirs[i] * norm(v3) / (4 * area);

      int pos = find_position(v, nodes);
      for (int k=0; k<d; k++) N(pos, k) += normal[k] * u; 
    }
  }

  for (int i=0; i<num_nodes; i++) {  // calculate R T R^T / volume
    for (int k=0; k<d; k++) v1[k] = N(i, k);
    v2 = T * v1;

    for (int j=i; j<num_nodes; j++) {
      for (int k=0; k<d; k++) v1[k] = N(j, k); 
      Ac(i, j) = (v1 * v2) / volume;
    }
  }

  const AmanziGeometry::Point& cm = mesh_->cell_centroid(cell);
  for (int i=0; i<num_nodes; i++) {
    int v = nodes[i];
    mesh_->node_get_coordinates(v, &p);
    for (int k=0; k<d; k++) N(i, k) = p[k] - cm[k];
    N(i, d) = 1;  // additional colum is added to the consistency condition
  }
  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* Consistency condition for stifness matrix in geomechanics. 
* Only the upper triangular part of Ac is calculated.
****************************************************************** */
int MFD3D::H1_consistency_elasticity(int cell,
                                     Tensor& T,
                                     Teuchos::SerialDenseMatrix<int, double>& N,
                                     Teuchos::SerialDenseMatrix<int, double>& Ac)
{
  return WHETSTONE_ELEMENTAL_MATRIX_OK;  // (lipnikov@lanl.gov)
}


/* ******************************************************************
* Conventional Gramm-Schimdt orthogonalization of colums of matrix N. 
****************************************************************** */
void MFD3D::stability_scalar(int cell,
                             Teuchos::SerialDenseMatrix<int, double>& N,
                             Teuchos::SerialDenseMatrix<int, double>& Mc,
                             Teuchos::SerialDenseMatrix<int, double>& M)
{
  gramm_schmidt(N);
 
  int nrows = Mc.numRows();
  int ncols = N.numCols();
  double volume = mesh_->cell_volume(cell);

  double scale = 0.0;
  for (int i=0; i<nrows; i++) scale += Mc(i, i);
  scale /= (nrows * volume);

  for (int i=0; i<nrows; i++) {
    for (int j=i; j<nrows; j++) M(i, j) = Mc(i, j);
  }

  for (int i=0; i<nrows; i++ ) {  // add projector (I - N^T N) to matrix M
    M(i, i) += scale;

    for (int j=i; j<nrows; j++) {
      double s = 0.0;
      for (int k=0; k<ncols; k++)  s += N(i, k) * N(j, k);
      M(i, j) -= s * scale;
    }
  }

  for (int i=0; i<nrows; i++) {  // symmetrization
    for (int j=i+1; j<nrows; j++) M(j, i) = M(i, j);
  }
}


/* ******************************************************************
* Conventional Gramm-Schimdt orthogonalization of colums of matrix N. 
****************************************************************** */
void MFD3D::gramm_schmidt(Teuchos::SerialDenseMatrix<int, double>& N)
{
  int nrows = N.numRows();
  int ncols = N.numCols();

  int i, j, k;
  for (i=0; i<ncols; i++ ) {
    double l22 = 0.0;
    for (k=0; k<nrows; k++) l22 += N(k, i) * N(k, i);

    l22 = 1.0 / sqrt(l22);
    for (k=0; k<nrows; k++) N(k, i) *= l22;

    for (j=i+1; j<ncols; j++) {
      double s = 0.0;
      for (k=0; k<nrows; k++) s += N(k, i) * N(k, j);
      for (k=0; k<nrows; k++) N(k, j) -= s * N(k, i);  // orthogonolize i and j
    }
  }
}


/* ******************************************************************
* Extension of Mesh API. 
****************************************************************** */
int MFD3D::cell_get_face_adj_cell(const int cell, const int face)
{
  AmanziMesh::Entity_ID_List cells; 
  mesh_->face_get_cells(face, AmanziMesh::USED, &cells);
  int ncells = cells.size();
 
  if (ncells == 2) {
    int c2 = cells[0];
    if (cell == c2) c2 = cells[1];
    return c2;
  }
  return -1;  
}


/* ******************************************************************
* Returns position of the number v in the list of nodes.  
****************************************************************** */
int MFD3D::find_position(int v, AmanziMesh::Entity_ID_List nodes)
{
  for (int i=0; i<nodes.size(); i++) {
    if (nodes[i] == v) return i;
  }
  return -1;
}

}  // namespace WhetStone
}  // namespace Amanzi



