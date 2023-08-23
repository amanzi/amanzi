/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Operators

  Base class for testing advection-diffusion problems with gravity:
    -div (K k (grad(p) - g) - v p) = f
  where g is the gravity vector pointing downward of axis z or
  axis y in two dimensions, K is diffusion tensor, v is velocity,
  and k is scalar correction to diffusion coefficient.

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
#include "MFD3D_Lagrange.hh"
#include "Mesh.hh"
#include "NumericalIntegration.hh"
#include "WhetStoneFunction.hh"

class AnalyticBase : public Amanzi::WhetStone::WhetStoneFunction {
 public:
  AnalyticBase(Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh)
    : mesh_(mesh), d_(mesh->getSpaceDimension()){};
  ~AnalyticBase(){};

  // analytic solution for diffusion problem with gravity
  // -- tensorial diffusion coefficient
  virtual Amanzi::WhetStone::Tensor
  TensorDiffusivity(const Amanzi::AmanziGeometry::Point& p, double t) = 0;

  // -- scalar diffusion coefficient
  virtual double ScalarDiffusivity(const Amanzi::AmanziGeometry::Point& p, double t) { return 1.0; }

  // -- analytic solution p
  virtual double pressure_exact(const Amanzi::AmanziGeometry::Point& p, double t) const = 0;

  // -- gradient of continuous velocity grad(h), where h = p - g z
  virtual Amanzi::AmanziGeometry::Point
  gradient_exact(const Amanzi::AmanziGeometry::Point& p, double t) = 0;

  // -- analytic advection function
  virtual Amanzi::AmanziGeometry::Point
  advection_exact(const Amanzi::AmanziGeometry::Point& p, double t) = 0;

  // -- source term
  virtual double source_exact(const Amanzi::AmanziGeometry::Point& p, double t) = 0;

  // derived quantity: Darcy velocity -K k grad(h)
  virtual Amanzi::AmanziGeometry::Point
  velocity_exact(const Amanzi::AmanziGeometry::Point& p, double t)
  {
    auto K = TensorDiffusivity(p, t);
    double kr = ScalarDiffusivity(p, t);
    auto grad = gradient_exact(p, t);
    return -(K * grad) * kr;
  }

  // interface to function
  virtual double Value(const Amanzi::AmanziGeometry::Point& p) const
  {
    return pressure_exact(p, 0.0);
  }

  // access
  int get_dimension() const { return d_; }

  // error calculation
  void
  ComputeCellError(Epetra_MultiVector& p, double t, double& pnorm, double& l2_err, double& inf_err);
  void
  ComputeFaceError(Epetra_MultiVector& u, double t, double& unorm, double& l2_err, double& inf_err);
  void ComputeNodeError(Epetra_MultiVector& p,
                        double t,
                        double& pnorm,
                        double& l2_err,
                        double& inf_err,
                        double& hnorm,
                        double& h1_err);
  void ComputeEdgeError(Epetra_MultiVector& p,
                        double t,
                        double& pnorm,
                        double& l2_err,
                        double& inf_err,
                        double& hnorm,
                        double& h1_err);
  void ComputeEdgeMomentsError(Epetra_MultiVector& p,
                               double t,
                               int ngauss,
                               double& pnorm,
                               double& l2_err,
                               double& inf_err);

  // utilities
  Amanzi::AmanziGeometry::Point face_normal_exterior(int f, bool* flag);

  // communications
  void GlobalOp(std::string op, double* val, int n);

 protected:
  Teuchos::RCP<const Amanzi::AmanziMesh::Mesh> mesh_;
  int d_;
};


/* ******************************************************************
* Error for cell-based fields
****************************************************************** */
inline void
AnalyticBase::ComputeCellError(Epetra_MultiVector& p,
                               double t,
                               double& pnorm,
                               double& l2_err,
                               double& inf_err)
{
  pnorm = 0.0;
  l2_err = 0.0;
  inf_err = 0.0;

  int ncells = mesh_->getNumEntities(Amanzi::AmanziMesh::Entity_kind::CELL,
                                     Amanzi::AmanziMesh::Parallel_type::OWNED);
  for (int c = 0; c < ncells; c++) {
    const Amanzi::AmanziGeometry::Point& xc = mesh_->getCellCentroid(c);
    double tmp = pressure_exact(xc, t);
    double volume = mesh_->getCellVolume(c);

    // std::cout << c << " xc=" << xc << " p: " << tmp << " " << p[0][c] << std::endl;
    l2_err += std::pow(tmp - p[0][c], 2.0) * volume;
    inf_err = std::max(inf_err, fabs(tmp - p[0][c]));
    pnorm += std::pow(tmp, 2.0) * volume;
  }
#ifdef HAVE_MPI
  GlobalOp("sum", &pnorm, 1);
  GlobalOp("sum", &l2_err, 1);
  GlobalOp("max", &inf_err, 1);
#endif
  pnorm = sqrt(pnorm);
  l2_err = sqrt(l2_err);
}


/* ******************************************************************
* Error for face-based fields
****************************************************************** */
inline void
AnalyticBase::ComputeFaceError(Epetra_MultiVector& u,
                               double t,
                               double& unorm,
                               double& l2_err,
                               double& inf_err)
{
  unorm = 0.0;
  l2_err = 0.0;
  inf_err = 0.0;

  int nfaces = mesh_->getNumEntities(Amanzi::AmanziMesh::Entity_kind::FACE,
                                     Amanzi::AmanziMesh::Parallel_type::OWNED);
  for (int f = 0; f < nfaces; f++) {
    double area = mesh_->getFaceArea(f);
    const Amanzi::AmanziGeometry::Point& normal = mesh_->getFaceNormal(f);
    const Amanzi::AmanziGeometry::Point& xf = mesh_->getFaceCentroid(f);
    const Amanzi::AmanziGeometry::Point& velocity = velocity_exact(xf, t);
    double tmp = velocity * normal;

    l2_err += std::pow((tmp - u[0][f]) / area, 2.0);
    inf_err = std::max(inf_err, fabs(tmp - u[0][f]) / area);
    unorm += std::pow(tmp / area, 2.0);
    // std::cout << f << " xf=" << xf << " u=" << u[0][f] << " u_ex=" << tmp << std::endl;
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
* Error for face-based fields
****************************************************************** */
inline void
AnalyticBase::ComputeNodeError(Epetra_MultiVector& p,
                               double t,
                               double& pnorm,
                               double& l2_err,
                               double& inf_err,
                               double& hnorm,
                               double& h1_err)
{
  pnorm = 0.0;
  l2_err = 0.0;
  inf_err = 0.0;
  hnorm = 0.0;
  h1_err = 0.0;

  Amanzi::AmanziGeometry::Point xv(d_);
  Amanzi::AmanziGeometry::Point grad(d_);

  Teuchos::ParameterList plist;
  plist.set<int>("method order", 1);
  Amanzi::WhetStone::MFD3D_Lagrange mfd(plist, mesh_);

  Amanzi::WhetStone::Polynomial poly(d_, 1);
  Amanzi::AmanziMesh::Entity_ID_List nodes;
  int ncells = mesh_->getNumEntities(Amanzi::AmanziMesh::Entity_kind::CELL,
                                     Amanzi::AmanziMesh::Parallel_type::OWNED);

  for (int c = 0; c < ncells; c++) {
    double volume = mesh_->getCellVolume(c);

    nodes = mesh_->getCellNodes(c);
    int nnodes = nodes.size();
    std::vector<Amanzi::WhetStone::Polynomial> cell_solution(nnodes);

    for (int k = 0; k < nnodes; k++) {
      int v = nodes[k];
      cell_solution[k].Reshape(d_, 0);
      cell_solution[k](0) = p[0][v];

      xv = mesh_->getNodeCoordinate(v);
      double tmp = pressure_exact(xv, t);

      if (std::abs(tmp - p[0][v]) > .01) {
        Amanzi::AmanziGeometry::Point xv2(2);
        xv2 = mesh_->getNodeCoordinate(v);
        // std::cout << v << " at " << xv << " error: " << tmp << " " << p[0][v] << std::endl;
      }
      l2_err += std::pow(tmp - p[0][v], 2.0) * volume / nnodes;
      inf_err = std::max(inf_err, fabs(tmp - p[0][v]));
      pnorm += std::pow(tmp, 2.0) * volume / nnodes;
    }

    const Amanzi::AmanziGeometry::Point& xc = mesh_->getCellCentroid(c);
    const Amanzi::AmanziGeometry::Point& grad_exact = gradient_exact(xc, t);
    mfd.L2Cell(c, cell_solution, cell_solution, NULL, poly);
    for (int k = 0; k < d_; ++k) grad[k] = poly(k + 1);

    h1_err += L22(grad - grad_exact) * volume;
    hnorm += L22(grad_exact) * volume;
  }
#ifdef HAVE_MPI
  GlobalOp("sum", &pnorm, 1);
  GlobalOp("sum", &l2_err, 1);
  GlobalOp("max", &inf_err, 1);
  GlobalOp("sum", &hnorm, 1);
  GlobalOp("sum", &h1_err, 1);
#endif
  pnorm = sqrt(pnorm);
  l2_err = sqrt(l2_err);

  hnorm = sqrt(hnorm);
  h1_err = sqrt(h1_err);
}


/* ******************************************************************
* Error for face-based fields
****************************************************************** */
inline void
AnalyticBase::ComputeEdgeError(Epetra_MultiVector& p,
                               double t,
                               double& pnorm,
                               double& l2_err,
                               double& inf_err,
                               double& hnorm,
                               double& h1_err)
{
  pnorm = 0.0;
  l2_err = 0.0;
  inf_err = 0.0;
  hnorm = 1.0; // missing code
  h1_err = 0.0;

  Amanzi::AmanziGeometry::Point grad(d_);

  Amanzi::AmanziMesh::Entity_ID_List edges;
  Amanzi::WhetStone::MFD3D_Diffusion mfd(mesh_);
  int ncells = mesh_->getNumEntities(Amanzi::AmanziMesh::Entity_kind::CELL,
                                     Amanzi::AmanziMesh::Parallel_type::OWNED);

  for (int c = 0; c < ncells; c++) {
    double volume = mesh_->getCellVolume(c);

    edges = mesh_->getCellEdges(c);
    int nedges = edges.size();
    std::vector<double> cell_solution(nedges);

    for (int k = 0; k < nedges; k++) {
      int e = edges[k];
      cell_solution[k] = p[0][e];

      const Amanzi::AmanziGeometry::Point& xe = mesh_->getEdgeCentroid(e);
      double tmp = pressure_exact(xe, t);
      l2_err += std::pow(tmp - p[0][e], 2.0) * volume / nedges;
      inf_err = std::max(inf_err, fabs(tmp - p[0][e]));
      pnorm += std::pow(tmp, 2.0) * volume / nedges;
      // std::cout << e << " at " << xe << " error: " << tmp << " " << p[0][e] << std::endl;
    }
  }
#ifdef HAVE_MPI
  GlobalOp("sum", &pnorm, 1);
  GlobalOp("sum", &l2_err, 1);
  GlobalOp("max", &inf_err, 1);
#endif
  pnorm = sqrt(pnorm);
  l2_err = sqrt(l2_err);
}


/* ******************************************************************
* Error in edge-based fields
****************************************************************** */
inline void
AnalyticBase::ComputeEdgeMomentsError(Epetra_MultiVector& p,
                                      double t,
                                      int ngauss,
                                      double& pnorm,
                                      double& l2_err,
                                      double& inf_err)
{
  pnorm = 0.0;
  l2_err = 0.0;
  inf_err = 0.0;

  Amanzi::AmanziMesh::Entity_ID n0, n1;
  Amanzi::AmanziGeometry::Point x0(d_), x1(d_), xv(d_);

  Amanzi::AmanziMesh::Entity_ID_List edges;
  Amanzi::WhetStone::MFD3D_Diffusion mfd(mesh_);
  int ncells = mesh_->getNumEntities(Amanzi::AmanziMesh::Entity_kind::CELL,
                                     Amanzi::AmanziMesh::Parallel_type::OWNED);

  for (int c = 0; c < ncells; c++) {
    double volume = mesh_->getCellVolume(c);

    edges = mesh_->getCellEdges(c);
    int nedges = edges.size();

    for (int k = 0; k < nedges; k++) {
      int e = edges[k];

      auto nodes = mesh_->getEdgeNodes(e);
      n0 = nodes[0];
      n1 = nodes[1];
      x0 = mesh_->getNodeCoordinate(n0);
      x1 = mesh_->getNodeCoordinate(n1);

      double s0(0.0);
      for (int n = 0; n <= ngauss; ++n) {
        double gp = Amanzi::WhetStone::q1d_points[ngauss - 1][n];
        double gw = Amanzi::WhetStone::q1d_weights[ngauss - 1][n];

        xv = x0 * gp + x1 * (1.0 - gp);
        s0 += gw * pressure_exact(xv, t);
        // s1 += gw * pressure_exact(xv, t) * (0.5 - gp);
      }

      l2_err += std::pow(s0 - p[0][e], 2.0) * volume / nedges;
      inf_err = std::max(inf_err, fabs(s0 - p[0][e]));
      pnorm += std::pow(s0, 2.0) * volume / nedges;
      // std::cout << e << " at " << (x0 + x1) / 2 << " error: " << s0 << " " << p[0][e] << std::endl;
    }
  }
#ifdef HAVE_MPI
  GlobalOp("sum", &pnorm, 1);
  GlobalOp("sum", &l2_err, 1);
  GlobalOp("max", &inf_err, 1);
#endif
  pnorm = sqrt(pnorm);
  l2_err = sqrt(l2_err);
}


/* ******************************************************************
* Exterior normal
****************************************************************** */
inline Amanzi::AmanziGeometry::Point
AnalyticBase::face_normal_exterior(int f, bool* flag)
{
  Amanzi::AmanziMesh::Entity_ID_List cells;
  cells = mesh_->getFaceCells(f, Amanzi::AmanziMesh::Parallel_type::ALL);
  *flag = (cells.size() == 1);

  int dir;
  Amanzi::AmanziGeometry::Point normal(d_);
  if (*flag) normal = mesh_->getFaceNormal(f, cells[0], &dir);

  return normal;
}


/* ******************************************************************
* Collective communications.
****************************************************************** */
inline void
AnalyticBase::GlobalOp(std::string op, double* val, int n)
{
  double* val_tmp = new double[n];
  for (int i = 0; i < n; ++i) val_tmp[i] = val[i];

  if (op == "sum")
    mesh_->getComm()->SumAll(val_tmp, val, n);
  else if (op == "max")
    mesh_->getComm()->MaxAll(val_tmp, val, n);

  delete[] val_tmp;
}

#endif
