/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Base class for testing diffusion problems with gravity:
    -div (K (grad(p) - g)) = f
  where g is the gravity vector pointing downward of axis z 
  or axis y in two dimensions.


  List of problems.  Note that all are 2D:

  Analytic00: polynomial solution with constant, scalar coefficient
  Analytic01: non-polynomial solution with full, non-constant tensor
  Analytic02: linear solution with constant, tensor coefficient
  Analytic03: non-polynomial solution with discontinuous (scalar)
              coefficient
  Analytic03b: same as 03, but with the coef as a scalar instead of
               scaling the tensor
  Analytic04: non-polynomial solution with non-constant scalar
              coefficient, coefficient can be zero
  Analytic05: linear solution with non-symmetric, non-constant
              tensor coefficient
*/

#ifndef AMANZI_OPERATOR_ANALYTIC_BASE_HH_
#define AMANZI_OPERATOR_ANALYTIC_BASE_HH_

#include "Epetra_MultiVector.h"

#include "MFD3D_Diffusion.hh"
#include "Mesh.hh"
#include "NumericalIntegration.hh"

class AnalyticBase {
 public:
  AnalyticBase(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh) : mesh_(mesh) {};
  ~AnalyticBase() {};

  // analytic solution for diffusion problem with gravity
  // -- diffusion tensor T
  virtual Amanzi::WhetStone::Tensor Tensor(const Amanzi::AmanziGeometry::Point& p, double t) = 0;

  // -- scalar component of the coefficient
  virtual double ScalarCoefficient(const Amanzi::AmanziGeometry::Point& p, double t) { return 1.0; }

  // -- analytic solution p
  virtual double pressure_exact(const Amanzi::AmanziGeometry::Point& p, double t) = 0;

  // -- gradient of continuous velocity grad(h), where h = p + g z
  virtual Amanzi::AmanziGeometry::Point gradient_exact(const Amanzi::AmanziGeometry::Point& p, double t) = 0;

  // -- source term
  virtual double source_exact(const Amanzi::AmanziGeometry::Point& p, double t) = 0;

  // derived quantity: Darcy velocity -K * grad(h)
  virtual Amanzi::AmanziGeometry::Point velocity_exact(const Amanzi::AmanziGeometry::Point& p, double t) {
    Amanzi::WhetStone::Tensor K = Tensor(p, t);
    Amanzi::AmanziGeometry::Point g = gradient_exact(p, t);
    double kr = ScalarCoefficient(p, t);
    return -(K * g) * kr;
  }

  // error calculation
  void ComputeCellError(Epetra_MultiVector& p, double t, double& pnorm, double& l2_err, double& inf_err) {
    pnorm = 0.0;
    l2_err = 0.0;
    inf_err = 0.0;

    int ncells = mesh_->num_entities(Amanzi::AmanziMesh::CELL, Amanzi::AmanziMesh::OWNED);
    for (int c = 0; c < ncells; c++) {
      const Amanzi::AmanziGeometry::Point& xc = mesh_->cell_centroid(c);
      double tmp = pressure_exact(xc, t);
      double volume = mesh_->cell_volume(c);

      // std::cout << c << " " << tmp << " " << p[0][c] << std::endl;
      l2_err += std::pow(tmp - p[0][c], 2.0) * volume;
      inf_err = std::max(inf_err, fabs(tmp - p[0][c]));
      pnorm += std::pow(tmp, 2.0) * volume;
    }
#ifdef HAVE_MPI
    double tmp = pnorm;
    mesh_->get_comm()->SumAll(&tmp, &pnorm, 1);
    tmp = l2_err;
    mesh_->get_comm()->SumAll(&tmp, &l2_err, 1);
    tmp = inf_err;
    mesh_->get_comm()->MaxAll(&tmp, &inf_err, 1);
#endif
    pnorm = sqrt(pnorm);
    l2_err = sqrt(l2_err);
  }

  void ComputeFaceError(Epetra_MultiVector& u, double t, double& unorm, double& l2_err, double& inf_err) {
    unorm = 0.0;
    l2_err = 0.0;
    inf_err = 0.0;

    int nfaces = mesh_->num_entities(Amanzi::AmanziMesh::FACE, Amanzi::AmanziMesh::OWNED);
    for (int f = 0; f < nfaces; f++) {
      double area = mesh_->face_area(f);
      const Amanzi::AmanziGeometry::Point& normal = mesh_->face_normal(f);
      const Amanzi::AmanziGeometry::Point& xf = mesh_->face_centroid(f);
      const Amanzi::AmanziGeometry::Point& velocity = velocity_exact(xf, t);
      double tmp = velocity * normal;

      l2_err += std::pow((tmp - u[0][f]) / area, 2.0);
      inf_err = std::max(inf_err, fabs(tmp - u[0][f]) / area);
      unorm += std::pow(tmp / area, 2.0);
      // std::cout << f << " " << xf << " " << u[0][f] << " " << tmp << std::endl;
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

  void ComputeNodeError(
      Epetra_MultiVector& p, double t,
      double& pnorm, double& l2_err, double& inf_err, double& hnorm, double& h1_err) {
    pnorm = 0.0;
    l2_err = 0.0;
    inf_err = 0.0;
    hnorm = 0.0;
    h1_err = 0.0;

    int d = mesh_->space_dimension();
    Amanzi::AmanziGeometry::Point xv(d);
    Amanzi::AmanziGeometry::Point grad(d);

    Amanzi::AmanziMesh::Entity_ID_List nodes;
    Amanzi::WhetStone::MFD3D_Diffusion mfd(mesh_);
    int ncells = mesh_->num_entities(Amanzi::AmanziMesh::CELL, Amanzi::AmanziMesh::OWNED);

    for (int c = 0; c < ncells; c++) {
      double volume = mesh_->cell_volume(c);

      mesh_->cell_get_nodes(c, &nodes);
      int nnodes = nodes.size();
      std::vector<double> cell_solution(nnodes);

      for (int k = 0; k < nnodes; k++) {
        int v = nodes[k];
        cell_solution[k] = p[0][v];

        mesh_->node_get_coordinates(v, &xv);
        double tmp = pressure_exact(xv, t);

        if (std::abs(tmp - p[0][v]) > .01) {
          Amanzi::AmanziGeometry::Point xv(2);
          mesh_->node_get_coordinates(v, &xv);
          // std::cout << v << " at " << xv << " error: " << tmp << " " << p[0][v] << std::endl;
        }
        l2_err += std::pow(tmp - p[0][v], 2.0) * volume / nnodes;
        inf_err = std::max(inf_err, fabs(tmp - p[0][v]));
        pnorm += std::pow(tmp, 2.0) * volume / nnodes;
      }

      const Amanzi::AmanziGeometry::Point& xc = mesh_->cell_centroid(c);
      const Amanzi::AmanziGeometry::Point& grad_exact = gradient_exact(xc, t);
      mfd.RecoverGradient_StiffnessMatrix(c, cell_solution, grad);

      h1_err += L22(grad - grad_exact) * volume;
      hnorm += L22(grad_exact) * volume;
    }
#ifdef HAVE_MPI
    double tmp = pnorm;
    mesh_->get_comm()->SumAll(&tmp, &pnorm, 1);
    tmp = l2_err;
    mesh_->get_comm()->SumAll(&tmp, &l2_err, 1);
    tmp = inf_err;
    mesh_->get_comm()->MaxAll(&tmp, &inf_err, 1);
    tmp = hnorm;
    mesh_->get_comm()->SumAll(&tmp, &hnorm, 1);
    tmp = h1_err;
    mesh_->get_comm()->SumAll(&tmp, &h1_err, 1);
#endif
    pnorm = sqrt(pnorm);
    l2_err = sqrt(l2_err);

    hnorm = sqrt(hnorm);
    h1_err = sqrt(h1_err);
  }

  void ComputeEdgeError(
      Epetra_MultiVector& p, double t,
      double& pnorm, double& l2_err, double& inf_err, double& hnorm, double& h1_err) {
    pnorm = 0.0;
    l2_err = 0.0;
    inf_err = 0.0;
    hnorm = 1.0;  // missing code
    h1_err = 0.0;

    int d = mesh_->space_dimension();
    Amanzi::AmanziGeometry::Point grad(d);

    Amanzi::AmanziMesh::Entity_ID_List edges;
    Amanzi::WhetStone::MFD3D_Diffusion mfd(mesh_);
    int ncells = mesh_->num_entities(Amanzi::AmanziMesh::CELL, Amanzi::AmanziMesh::OWNED);

    for (int c = 0; c < ncells; c++) {
      double volume = mesh_->cell_volume(c);

      mesh_->cell_get_edges(c, &edges);
      int nedges = edges.size();
      std::vector<double> cell_solution(nedges);

      for (int k = 0; k < nedges; k++) {
        int e = edges[k];
        cell_solution[k] = p[0][e];

        const Amanzi::AmanziGeometry::Point& xe = mesh_->edge_centroid(e);
        double tmp = pressure_exact(xe, t);
        l2_err += std::pow(tmp - p[0][e], 2.0) * volume / nedges;
        inf_err = std::max(inf_err, fabs(tmp - p[0][e]));
        pnorm += std::pow(tmp, 2.0) * volume / nedges;
        // std::cout << e << " at " << xe << " error: " << tmp << " " << p[0][e] << std::endl;
      }
    }
#ifdef HAVE_MPI
    double tmp = pnorm;
    mesh_->get_comm()->SumAll(&tmp, &pnorm, 1);
    tmp = l2_err;
    mesh_->get_comm()->SumAll(&tmp, &l2_err, 1);
    tmp = inf_err;
    mesh_->get_comm()->MaxAll(&tmp, &inf_err, 1);
#endif
    pnorm = sqrt(pnorm);
    l2_err = sqrt(l2_err);
  }

  void ComputeEdgeMomentsError(
      Epetra_MultiVector& p, double t, int ngauss, 
      double& pnorm, double& l2_err, double& inf_err) {
    pnorm = 0.0;
    l2_err = 0.0;
    inf_err = 0.0;

    int d = mesh_->space_dimension();

    Amanzi::AmanziMesh::Entity_ID n0, n1;
    Amanzi::AmanziGeometry::Point x0(d), x1(d), xv(d);

    Amanzi::AmanziMesh::Entity_ID_List edges;
    Amanzi::WhetStone::MFD3D_Diffusion mfd(mesh_);
    int ncells = mesh_->num_entities(Amanzi::AmanziMesh::CELL, Amanzi::AmanziMesh::OWNED);

    for (int c = 0; c < ncells; c++) {
      double volume = mesh_->cell_volume(c);

      mesh_->cell_get_edges(c, &edges);
      int nedges = edges.size();

      for (int k = 0; k < nedges; k++) {
        int e = edges[k];

        mesh_->edge_get_nodes(e, &n0, &n1);
        mesh_->node_get_coordinates(n0, &x0);
        mesh_->node_get_coordinates(n1, &x1);

        double s0(0.0), s1(0.0);
        for (int n = 0; n <= ngauss; ++n) {
          double gp = Amanzi::WhetStone::q1d_points[ngauss - 1][n];
          double gw = Amanzi::WhetStone::q1d_weights[ngauss - 1][n];

          xv = x0 * gp + x1 * (1.0 - gp);
          s0 += gw * pressure_exact(xv, t);
          s1 += gw * pressure_exact(xv, t) * (0.5 - gp);
        } 

        l2_err += std::pow(s0 - p[0][e], 2.0) * volume / nedges;
        inf_err = std::max(inf_err, fabs(s0 - p[0][e]));
        pnorm += std::pow(s0, 2.0) * volume / nedges;
        // std::cout << e << " at " << (x0 + x1) / 2 << " error: " << s0 << " " << p[0][e] << std::endl;
      }
    }
#ifdef HAVE_MPI
    double tmp = pnorm;
    mesh_->get_comm()->SumAll(&tmp, &pnorm, 1);
    tmp = l2_err;
    mesh_->get_comm()->SumAll(&tmp, &l2_err, 1);
    tmp = inf_err;
    mesh_->get_comm()->MaxAll(&tmp, &inf_err, 1);
#endif
    pnorm = sqrt(pnorm);
    l2_err = sqrt(l2_err);
  }

 protected:
  Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh_;
};

#endif

