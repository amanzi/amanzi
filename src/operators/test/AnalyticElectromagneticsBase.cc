/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Operators

  Base class for verifying numerical schemes for Maxwell's equations.
*/

#include "Epetra_MultiVector.h"

#include "AnalyticElectromagneticsBase.hh"

/* ******************************************************************
* Error calculation of faces.
****************************************************************** */
void
AnalyticElectromagneticsBase::ComputeFaceError(Epetra_MultiVector& u,
                                               double t,
                                               double& unorm,
                                               double& l2_err,
                                               double& inf_err)
{
  unorm = 0.0;
  l2_err = 0.0;
  inf_err = 0.0;

  int nfaces =
    mesh_->getNumEntities(Amanzi::AmanziMesh::Entity_kind::FACE, Amanzi::AmanziMesh::Parallel_type::OWNED);
  for (int f = 0; f < nfaces; f++) {
    double area = mesh_->getFaceArea(f);
    const Amanzi::AmanziGeometry::Point& normal = mesh_->getFaceNormal(f);
    const Amanzi::AmanziGeometry::Point& xf = mesh_->getFaceCentroid(f);
    const Amanzi::AmanziGeometry::Point& B = magnetic_exact(xf, t);
    double tmp = (normal * B) / area;

    l2_err += std::pow((tmp - u[0][f]), 2.0);
    inf_err = std::max(inf_err, fabs(tmp - u[0][f]));
    unorm += std::pow(tmp, 2.0);
    // std::cout << f << " x=" << xf << "  normal=" << normal << " " << u[0][f] << " " << tmp << std::endl;
  }
#ifdef HAVE_MPI
  double tmp = unorm;
  mesh_->getComm()->SumAll(&tmp, &unorm, 1);
  tmp = l2_err;
  mesh_->getComm()->SumAll(&tmp, &l2_err, 1);
  tmp = inf_err;
  mesh_->getComm()->MaxAll(&tmp, &inf_err, 1);
#endif
  unorm = sqrt(unorm);
  l2_err = sqrt(l2_err);
};


/* ******************************************************************
* Error calculation of edges.
****************************************************************** */
void
AnalyticElectromagneticsBase::ComputeEdgeError(Epetra_MultiVector& u,
                                               double t,
                                               double& unorm,
                                               double& l2_err,
                                               double& inf_err)
{
  unorm = 0.0;
  l2_err = 0.0;
  inf_err = 0.0;

  int nedges =
    mesh_->getNumEntities(Amanzi::AmanziMesh::Entity_kind::EDGE, Amanzi::AmanziMesh::Parallel_type::OWNED);
  for (int e = 0; e < nedges; e++) {
    double len = mesh_->getEdgeLength(e);
    const Amanzi::AmanziGeometry::Point& tau = mesh_->getEdgeVector(e);
    const Amanzi::AmanziGeometry::Point& xe = mesh_->getEdgeCentroid(e);

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
void
AnalyticElectromagneticsBase::ComputeNodeError(Epetra_MultiVector& u,
                                               double t,
                                               double& unorm,
                                               double& l2_err,
                                               double& inf_err)
{
  unorm = 0.0;
  l2_err = 0.0;
  inf_err = 0.0;

  int d = mesh_->getSpaceDimension();
  Amanzi::AmanziGeometry::Point xv(d);

  Amanzi::AmanziMesh::Entity_ID_List nodes;
  int ncells =
    mesh_->getNumEntities(Amanzi::AmanziMesh::Entity_kind::CELL, Amanzi::AmanziMesh::Parallel_type::OWNED);

  for (int c = 0; c < ncells; c++) {
    double volume = mesh_->getCellVolume(c);

    nodes = mesh_->getCellNodes(c);
    int nnodes = nodes.size();

    for (int k = 0; k < nnodes; k++) {
      int v = nodes[k];
      xv = mesh_->getNodeCoordinate(v);
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
void
AnalyticElectromagneticsBase::GlobalOp(std::string op, double* val, int n)
{
  double* val_tmp = new double[n];
  for (int i = 0; i < n; ++i) val_tmp[i] = val[i];

  if (op == "sum")
    mesh_->getComm()->SumAll(val_tmp, val, n);
  else if (op == "max")
    mesh_->getComm()->MaxAll(val_tmp, val, n);

  delete[] val_tmp;
}
