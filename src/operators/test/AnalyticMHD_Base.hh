/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Base class for testing MHD problems.
*/

#ifndef AMANZI_OPERATOR_ANALYTIC_MHD_BASE_HH_
#define AMANZI_OPERATOR_ANALYTIC_MHD_BASE_HH_

#include "mfd3d_electromagnetics.hh"
#include "Mesh.hh"

class AnalyticMHD_Base {
 public:
  AnalyticMHD_Base(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh) : mesh_(mesh) {};
  ~AnalyticMHD_Base() {};

  // analytic solution for MHD problem
  // -- resitivity tensor T
  virtual Amanzi::WhetStone::Tensor Tensor(const Amanzi::AmanziGeometry::Point& p, double t) = 0;
  // -- analytic solution E
  virtual Amanzi::AmanziGeometry::Point electric_exact(const Amanzi::AmanziGeometry::Point& p, double t) = 0;
  virtual Amanzi::AmanziGeometry::Point magnetic_exact(const Amanzi::AmanziGeometry::Point& p, double t) = 0;
  // -- source term
  virtual Amanzi::AmanziGeometry::Point source_exact(const Amanzi::AmanziGeometry::Point& p, double t) = 0;

  // error calculation
  void ComputeEdgeError(Epetra_MultiVector& u, double t, double& unorm, double& l2_err, double& inf_err) {
    unorm = 0.0;
    l2_err = 0.0;
    inf_err = 0.0;

    int n1, n2; 
    Amanzi::AmanziGeometry::Point p1(3), p2(3), xe(3);

    int nedges = mesh_->num_entities(Amanzi::AmanziMesh::EDGE, Amanzi::AmanziMesh::OWNED);
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
      // std::cout << e << " " << u[0][e] << " " << tmp << std::endl;
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
  }

 protected:
  Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh_;
};

#endif

