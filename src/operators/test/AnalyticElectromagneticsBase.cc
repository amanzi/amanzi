/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Base class for verifying numerical schemes for Maxwell's equations.
*/

#include "Epetra_MultiVector.h"

#include "AnalyticElectromagneticsBase.hh"

/* ******************************************************************
* Error calculation of faces.
****************************************************************** */
void AnalyticElectromagneticsBase::ComputeFaceError(
    Epetra_MultiVector& u, double t, double& unorm, double& l2_err, double& inf_err)
{
  unorm = 0.0;
  l2_err = 0.0;
  inf_err = 0.0;

  int nfaces = mesh_->num_entities(Amanzi::AmanziMesh::FACE, Amanzi::AmanziMesh::Parallel_type::OWNED);
  for (int f = 0; f < nfaces; f++) {
    double area = mesh_->face_area(f);
    const Amanzi::AmanziGeometry::Point& normal = mesh_->face_normal(f);
    const Amanzi::AmanziGeometry::Point& xf = mesh_->face_centroid(f);
    const Amanzi::AmanziGeometry::Point& B = magnetic_exact(xf, t);
    double tmp = (normal * B) / area;

    l2_err += std::pow((tmp - u[0][f]), 2.0);
    inf_err = std::max(inf_err, fabs(tmp - u[0][f]));
    unorm += std::pow(tmp, 2.0);
    // std::cout << f << " x=" << xf << "  normal=" << normal << " " << u[0][f] << " " << tmp << std::endl;
  }
#ifdef HAVE_MPI
  double tmp = unorm;
  mesh_->get_comm()->SumAll(&tmp, &unorm, 1);
  tmp = l2_err;
  mesh_->get_comm()->SumAll(&tmp, &l2_err, 1);
  tmp = inf_err;
  mesh_->get_comm()->MaxAll(&tmp, &inf_err, 1);
#endif
  unorm = sqrt(unorm);
  l2_err = sqrt(l2_err);
};


/* ******************************************************************
* Error calculation of edges.
****************************************************************** */
void AnalyticElectromagneticsBase::ComputeEdgeError(
    Epetra_MultiVector& u, double t, double& unorm, double& l2_err, double& inf_err)
{
  unorm = 0.0;
  l2_err = 0.0;
  inf_err = 0.0;

  int nedges = mesh_->num_entities(Amanzi::AmanziMesh::EDGE, Amanzi::AmanziMesh::Parallel_type::OWNED);
  for (int e = 0; e < nedges; e++) {
    double len = mesh_->edge_length(e);
    const Amanzi::AmanziGeometry::Point& tau = mesh_->edge_vector(e);
    const Amanzi::AmanziGeometry::Point& xe = mesh_->edge_centroid(e);

    const Amanzi::AmanziGeometry::Point& E = electric_exact(xe, t);
    double tmp = (E * tau) / len;

    l2_err += std::pow(tmp - u[0][e], 2.0);
    inf_err = std::max(inf_err, fabs(tmp - u[0][e]));
    unorm += std::pow(tmp, 2.0);
    // std::cout << e << " xe=" << xe << " E=" << u[0][e] << " Eex=" << tmp << " L2err=" << l2_err << std::endl;
  }
#ifdef HAVE_MPI
  GlobalOp("sum", &unorm, 1);
  GlobalOp("sum", &l2_err, 1);
  GlobalOp("max", &inf_err, 1);
#endif
  unorm = sqrt(unorm);
  l2_err = sqrt(l2_err);
};


/* ******************************************************************
* Error calculation of nodes.
****************************************************************** */
void AnalyticElectromagneticsBase::ComputeNodeError(
    Epetra_MultiVector& u, double t, double& unorm, double& l2_err, double& inf_err)
{
  unorm = 0.0;
  l2_err = 0.0;
  inf_err = 0.0;

  int d = mesh_->space_dimension();
  Amanzi::AmanziGeometry::Point xv(d);

  Amanzi::AmanziMesh::Entity_ID_List nodes;
  int ncells = mesh_->num_entities(Amanzi::AmanziMesh::CELL, Amanzi::AmanziMesh::Parallel_type::OWNED);

  for (int c = 0; c < ncells; c++) {
    double volume = mesh_->cell_volume(c);

    mesh_->cell_get_nodes(c, &nodes);
    int nnodes = nodes.size();

    for (int k = 0; k < nnodes; k++) {
      int v = nodes[k];
      mesh_->node_get_coordinates(v, &xv);
      double tmp = (electric_exact(xv, t))[2];

      // std::cout << v << " at " << xv << " error: " << tmp << " " << u[0][v] << std::endl;
      l2_err += std::pow(tmp - u[0][v], 2.0) * volume / nnodes;
      inf_err = std::max(inf_err, fabs(tmp - u[0][v]));
      unorm += std::pow(tmp, 2.0) * volume / nnodes;
    }
  }
#ifdef HAVE_MPI
  GlobalOp("sum", &unorm, 1);
  GlobalOp("sum", &l2_err, 1);
  GlobalOp("max", &inf_err, 1);
#endif
  unorm = sqrt(unorm);
  l2_err = sqrt(l2_err);
}


/* ******************************************************************
* Collective communications.
****************************************************************** */
void AnalyticElectromagneticsBase::GlobalOp(std::string op, double* val, int n)
{
  double* val_tmp = new double[n];
  for (int i = 0; i < n; ++i) val_tmp[i] = val[i];

  if (op == "sum") 
    mesh_->get_comm()->SumAll(val_tmp, val, n);
  else if (op == "max") 
    mesh_->get_comm()->MaxAll(val_tmp, val, n);

  delete[] val_tmp;
}
