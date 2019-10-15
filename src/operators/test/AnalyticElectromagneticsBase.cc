/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Konstantin Lipnikov (lipnikov@lanl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>

#include "Epetra_MultiVector.h"

#include "AnalyticElectromagneticsBase.hh"

/* ******************************************************************
 * Error calculation of faces.
 ****************************************************************** */
void
AnalyticElectromagneticsBase::ComputeFaceError(Epetra_MultiVector& u, double t,
                                               double& unorm, double& l2_err,
                                               double& inf_err)
{
  unorm = 0.0;
  l2_err = 0.0;
  inf_err = 0.0;

  int nfaces = mesh_->num_entities(Amanzi::AmanziMesh::FACE,
                                   Amanzi::AmanziMesh::Parallel_type::OWNED);
  for (int f = 0; f < nfaces; f++) {
    double area = mesh_->face_area(f);
    const Amanzi::AmanziGeometry::Point& normal = mesh_->face_normal(f);
    const Amanzi::AmanziGeometry::Point& xf = mesh_->face_centroid(f);
    const Amanzi::AmanziGeometry::Point& B = magnetic_exact(xf, t);
    double tmp = (normal * B) / area;

    l2_err += std::pow((tmp - u[0][f]), 2.0);
    inf_err = std::max(inf_err, fabs(tmp - u[0][f]));
    unorm += std::pow(tmp, 2.0);
    // std::cout << f << " x=" << xf << "  normal=" << normal << " " << u[0][f]
    // << " " << tmp << std::endl;
  }
#ifdef HAVE_MPI
  double tmp = unorm;
  Teuchos::reduceAll(*mesh_->get_comm(), Teuchos::REDUCE_SUM, 1, &tmp, &unorm);
  tmp = l2_err;
  Teuchos::reduceAll(*mesh_->get_comm(), Teuchos::REDUCE_SUM, 1, &tmp, &l2_err);
  tmp = inf_err;
  Teuchos::reduceAll(
    *mesh_->get_comm(), Teuchos::REDUCE_MAX, 1, &tmp, &inf_err);
#endif
  unorm = sqrt(unorm);
  l2_err = sqrt(l2_err);
};


/* ******************************************************************
 * Error calculation of edges.
 ****************************************************************** */
void
AnalyticElectromagneticsBase::ComputeEdgeError(Epetra_MultiVector& u, double t,
                                               double& unorm, double& l2_err,
                                               double& inf_err)
{
  unorm = 0.0;
  l2_err = 0.0;
  inf_err = 0.0;

  int n1, n2;
  Amanzi::AmanziGeometry::Point p1(3), p2(3), xe(3);

  int nedges = mesh_->num_entities(Amanzi::AmanziMesh::EDGE,
                                   Amanzi::AmanziMesh::Parallel_type::OWNED);
  for (int e = 0; e < nedges; e++) {
    double len = mesh_->edge_length(e);
    const Amanzi::AmanziGeometry::Point& tau = mesh_->edge_vector(e);

    mesh_->edge_get_nodes(e, &n1, &n2);
    mesh_->node_get_coordinates(n1, &p1);
    mesh_->node_get_coordinates(n2, &p2);
    xe = (p1 + p2) / 2;

    const Amanzi::AmanziGeometry::Point& E = electric_exact(xe, t);
    double tmp = (E * tau) / len;

    l2_err += std::pow(tmp - u[0][e], 2.0);
    inf_err = std::max(inf_err, fabs(tmp - u[0][e]));
    unorm += std::pow(tmp, 2.0);
    // std::cout << e << " xe=" << xe << " E=" << u[0][e] << " Eex=" << tmp << "
    // L2err=" << l2_err << std::endl;
  }
#ifdef HAVE_MPI
  double tmp = unorm;
  Teuchos::reduceAll(*mesh_->get_comm(), Teuchos::REDUCE_SUM, 1, &tmp, &unorm);
  tmp = l2_err;
  Teuchos::reduceAll(*mesh_->get_comm(), Teuchos::REDUCE_SUM, 1, &tmp, &l2_err);
  tmp = inf_err;
  Teuchos::reduceAll(
    *mesh_->get_comm(), Teuchos::REDUCE_MAX, 1, &tmp, &inf_err);
#endif
  unorm = sqrt(unorm);
  l2_err = sqrt(l2_err);
};


/* ******************************************************************
 * Error calculation of nodes.
 ****************************************************************** */
void
AnalyticElectromagneticsBase::ComputeNodeError(Epetra_MultiVector& u, double t,
                                               double& unorm, double& l2_err,
                                               double& inf_err)
{
  unorm = 0.0;
  l2_err = 0.0;
  inf_err = 0.0;

  int d = mesh_->space_dimension();
  Amanzi::AmanziGeometry::Point xv(d);

  Amanzi::AmanziMesh::Entity_ID_List nodes;
  int ncells = mesh_->num_entities(Amanzi::AmanziMesh::CELL,
                                   Amanzi::AmanziMesh::Parallel_type::OWNED);

  for (int c = 0; c < ncells; c++) {
    double volume = mesh_->cell_volume(c, false);

    mesh_->cell_get_nodes(c, &nodes);
    int nnodes = nodes.size();

    for (int k = 0; k < nnodes; k++) {
      int v = nodes[k];
      mesh_->node_get_coordinates(v, &xv);
      double tmp = (electric_exact(xv, t))[2];

      // std::cout << v << " at " << xv << " error: " << tmp << " " << u[0][v]
      // << std::endl;
      l2_err += std::pow(tmp - u[0][v], 2.0) * volume / nnodes;
      inf_err = std::max(inf_err, fabs(tmp - u[0][v]));
      unorm += std::pow(tmp, 2.0) * volume / nnodes;
    }
  }
#ifdef HAVE_MPI
  double tmp = unorm;
  Teuchos::reduceAll(*mesh_->get_comm(), Teuchos::REDUCE_SUM, 1, &tmp, &unorm);
  tmp = l2_err;
  Teuchos::reduceAll(*mesh_->get_comm(), Teuchos::REDUCE_SUM, 1, &tmp, &l2_err);
  tmp = inf_err;
  Teuchos::reduceAll(
    *mesh_->get_comm(), Teuchos::REDUCE_MAX, 1, &tmp, &inf_err);
#endif
  unorm = sqrt(unorm);
  l2_err = sqrt(l2_err);
}
