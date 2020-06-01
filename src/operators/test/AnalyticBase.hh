/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Base class for testing advection-diffusion problems with gravity:
    -div (K k (grad(p) - g) - v p) = f
  where g is the gravity vector pointing downward of axis z or
  axis y in two dimensions, K is diffusion tensor, v is velocity,
  and k is scalar correction to diffusion coefficient.

  List of problems.  Note that all are 2D:

  Analytic00: polynomial solution with constant, scalar coefficient
  Analytic00b: 3D variant of Analytic00
  Analytic00c: Analytic00b with gravity
  Analytic01: non-polynomial solution with full, non-constant tensor
  Analytic01b: 3D variant of Analytic01
  Analytic02: linear solution with constant, tensor coefficient
  Analytic03: non-polynomial solution with discontinuous (scalar)
              coefficient
  Analytic03b: same as 03, but with the coef as a scalar instead of
               scaling the tensor
  Analytic04: non-polynomial solution with non-constant scalar
              coefficient, coefficient can be zero
  Analytic05: linear solution with non-symmetric, non-constant
              tensor coefficient
  Analytic06: ??? (Description missing)
*/

//! <MISSING_ONELINE_DOCSTRING>

#ifndef AMANZI_OPERATOR_ANALYTIC_BASE_HH_
#define AMANZI_OPERATOR_ANALYTIC_BASE_HH_

#include "Teuchos_CommHelpers.hpp"

#include "AmanziTypes.hh"
#include "AmanziVector.hh"
#include "CompositeVector.hh"
#include "MFD3D_Diffusion.hh"
#include "MFD3D_Lagrange.hh"
#include "Mesh.hh"
//#include "NumericalIntegration.hh"
#include "Quadrature1D.hh"
#include "Quadrature2D.hh"
#include "Quadrature3D.hh"

using namespace Amanzi;

class AnalyticBase { //: public WhetStone::WhetStoneFunction {
 public:
  AnalyticBase(int d) 
    : d_(d) {};
  virtual ~AnalyticBase() {};

  virtual std::string name() const = 0;
  int dimension() const { return d_; }

  // analytic solution for diffusion problem with gravity
  // -- tensorial diffusion coefficient
  virtual WhetStone::Tensor<Kokkos::HostSpace>
  TensorDiffusivity(const AmanziGeometry::Point& p, double t) const {
    WhetStone::Tensor<Kokkos::HostSpace> K(d_, 1);
    K(0,0) = 1.0;
    return K;
  }
    

  // -- scalar diffusion coefficient
  virtual double
  ScalarDiffusivity(const AmanziGeometry::Point& p, double t) const {
    return 1.0;
  }

  // -- analytic solution p
  virtual double pressure_exact(const AmanziGeometry::Point& p, double t) const = 0;
  double pressure_exact(const AmanziGeometry::Point& p) const {
    return pressure_exact(p, 0.0);
  }

  // -- gradient of continuous velocity grad(h), where h = p - g z
  virtual AmanziGeometry::Point
  gradient_exact(const AmanziGeometry::Point& p, double t) const = 0;

  // -- analytic advection function
  virtual AmanziGeometry::Point
  advection_exact(const AmanziGeometry::Point& p, double t) const = 0;

  // -- source term
  virtual double
  source_exact(const AmanziGeometry::Point& p, double t) const = 0;

  // derived quantity: Darcy velocity -K k grad(h)
  virtual AmanziGeometry::Point
  velocity_exact(const AmanziGeometry::Point& p, double t) const {
    auto K = TensorDiffusivity(p, t);
    double kr = ScalarDiffusivity(p, t);
    auto grad = gradient_exact(p, t);
    return -(K * grad) * kr;
  }

 protected:
  int d_;
};

//
// Nonmember helper functions
//

/* ******************************************************************
* Exterior normal
****************************************************************** */
inline
AmanziGeometry::Point
FaceNormalExterior(const AmanziMesh::Mesh& mesh, int f, bool* flag)
{
  AmanziMesh::Entity_ID_View cells;
  mesh.face_get_cells(f, AmanziMesh::Parallel_type::ALL, cells);
  *flag = (cells.extent(0) == 1);

  int dir;
  AmanziGeometry::Point normal;
  if (*flag) 
    normal = mesh.face_normal(f, false, cells[0], &dir);

  return normal;
}


/* ******************************************************************
* Collective communications.
****************************************************************** */
inline
void GlobalOp(const Comm_type& comm, std::string op, double* val, int n)
{
  double* val_tmp = new double[n];
  memcpy(val_tmp,val,n*sizeof(double)); 

  if (op == "sum")
    Teuchos::reduceAll(comm,Teuchos::REDUCE_SUM, 
      *val_tmp, Teuchos::outArg(*val));
  else if (op == "max")
    Teuchos::reduceAll(comm,Teuchos::REDUCE_MAX, 
      *val_tmp, Teuchos::outArg(*val));

  delete[] val_tmp;
}


/* ******************************************************************
* Error for cell-based fields
****************************************************************** */
inline
void ComputeCellError(const AnalyticBase& ana,
                      const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                      const CompositeVector& p_vec,
                      double t,
                      double& pnorm, double& l2_err, double& inf_err)
{
  pnorm = 0.0;
  l2_err = 0.0;
  inf_err = 0.0;

  int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);
  auto p = p_vec.ViewComponent<MirrorHost>("cell", false);
  for (int c = 0; c < ncells; c++) {
    const AmanziGeometry::Point& xc = mesh->cell_centroid(c);
    double tmp = ana.pressure_exact(xc, t);
    double volume = mesh->cell_volume(c);

    // std::cout << c << " xc=" << xc << " p: " << tmp << " " << p(c,0) <<
    // std::endl;
    l2_err += std::pow(tmp - p(c,0), 2.0) * volume;
    inf_err = std::max(inf_err, fabs(tmp - p(c,0)));
    pnorm += std::pow(tmp, 2.0) * volume;
    if (tmp - p(c,0) > 1.e-2)
      std::cout << c << " xc=" << xc << " p=" << p(c,0) << " p_ex=" << tmp << std::endl;
  }
#ifdef HAVE_MPI
  GlobalOp(*mesh->get_comm(), "sum", &pnorm, 1);
  GlobalOp(*mesh->get_comm(), "sum", &l2_err, 1);
  GlobalOp(*mesh->get_comm(), "max", &inf_err, 1);
#endif
  pnorm = sqrt(pnorm);
  l2_err = sqrt(l2_err);
}


/* ******************************************************************
* Error for face-based fields
****************************************************************** */
inline
void ComputeFaceError(const AnalyticBase& ana,
        const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
        const CompositeVector& u_vec, double t,
        double& unorm, double& l2_err, double& inf_err)
{
  unorm = 0.0;
  l2_err = 0.0;
  inf_err = 0.0;

  auto u = u_vec.ViewComponent<MirrorHost>("face", false);
  int nfaces = mesh->num_entities(AmanziMesh::FACE, AmanziMesh::Parallel_type::OWNED);
  for (int f = 0; f < nfaces; f++) {
    double area = mesh->face_area(f);
    const AmanziGeometry::Point& normal = mesh->face_normal(f);
    const AmanziGeometry::Point& xf = mesh->face_centroid(f);
    const AmanziGeometry::Point& velocity = ana.velocity_exact(xf, t);
    double tmp = velocity * normal;

    l2_err += std::pow((tmp - u(f,0)) / area, 2.0);
    inf_err = std::max(inf_err, fabs(tmp - u(f,0)) / area);
    unorm += std::pow(tmp / area, 2.0);
    if ((tmp - u(f,0)) / area > 1.e-1)
      std::cout << f << " xf=" << xf << " q=" << u(f,0) << " q_ex=" << tmp << " velocity=" << velocity << std::endl;
  }
#ifdef HAVE_MPI
  GlobalOp(*mesh->get_comm(), "sum", &unorm, 1);
  GlobalOp(*mesh->get_comm(), "sum", &l2_err, 1);
  GlobalOp(*mesh->get_comm(), "max", &inf_err, 1);
#endif
  unorm = sqrt(unorm);
  l2_err = sqrt(l2_err);
}


/* ******************************************************************
* Error for face-based fields
****************************************************************** */
inline
void ComputeNodeError(const AnalyticBase& ana,
        const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
        const CompositeVector& p_vec, double t,
        double& pnorm, double& l2_err, double& inf_err, double& hnorm, double& h1_err)
{
  pnorm = 0.0;
  l2_err = 0.0;
  inf_err = 0.0;
  hnorm = 0.0;
  h1_err = 0.0;

  AmanziGeometry::Point xv(ana.dimension());
  AmanziGeometry::Point grad(ana.dimension());

  Teuchos::ParameterList plist;
  plist.set<int>("method order", 1);
  WhetStone::MFD3D_Lagrange mfd(plist, mesh);

  WhetStone::Polynomial<> poly(ana.dimension(), 1);
  AmanziMesh::Entity_ID_List nodes;
  int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);

  auto p = p_vec.ViewComponent<MirrorHost>("node", false);
  for (int c = 0; c < ncells; c++) {
    double volume = mesh->cell_volume(c);

    mesh->cell_get_nodes(c, nodes);
    int nnodes = nodes.size();
    std::vector<WhetStone::Polynomial<>> cell_solution(nnodes);
    std::vector<WhetStone::Polynomial<>> vf_solutions(nnodes);

    for (int k = 0; k < nnodes; k++) {
      int v = nodes[k];
      cell_solution[k].reshape(ana.dimension(), 0);
      cell_solution[k](0) = p(v,0);

      mesh->node_get_coordinates(v, &xv);
      double tmp = ana.pressure_exact(xv, t);

      if (std::abs(tmp - p(v,0)) > .01) {
        AmanziGeometry::Point xv2(2);
        mesh->node_get_coordinates(v, &xv2);
        // std::cout << v << " at " << xv << " error: " << tmp << " " << p(v,0) << std::endl;
      }
      l2_err += std::pow(tmp - p(v,0), 2.0) * volume / nnodes;
      inf_err = std::max(inf_err, fabs(tmp - p(v,0)));
      pnorm += std::pow(tmp, 2.0) * volume / nnodes;
    }

    const AmanziGeometry::Point& xc = mesh->cell_centroid(c);
    const AmanziGeometry::Point& grad_exact = ana.gradient_exact(xc, t);
    mfd.L2Cell(c, cell_solution, vf_solutions, NULL, poly);
    for (int k = 0; k < ana.dimension(); ++k) grad[k] = poly(k + 1);

    h1_err += L22(grad - grad_exact) * volume;
    hnorm += L22(grad_exact) * volume;
  }
#ifdef HAVE_MPI
  GlobalOp(*mesh->get_comm(), "sum", &pnorm, 1);
  GlobalOp(*mesh->get_comm(), "sum", &l2_err, 1);
  GlobalOp(*mesh->get_comm(), "max", &inf_err, 1);
  GlobalOp(*mesh->get_comm(), "sum", &hnorm, 1);
  GlobalOp(*mesh->get_comm(), "sum", &h1_err, 1);
#endif
  pnorm = sqrt(pnorm);
  l2_err = sqrt(l2_err);

  hnorm = sqrt(hnorm);
  h1_err = sqrt(h1_err);
}


/* ******************************************************************
* Error for face-based fields
****************************************************************** */
inline
void ComputeEdgeError(const AnalyticBase& ana,
        const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,    
                      const CompositeVector& p_vec, double t,
                      double& pnorm, double& l2_err, double& inf_err,
                      double& hnorm, double& h1_err)
{
  pnorm = 0.0;
  l2_err = 0.0;
  inf_err = 0.0;
  hnorm = 1.0; // missing code
  h1_err = 0.0;

  AmanziGeometry::Point grad(ana.dimension());

  AmanziMesh::Entity_ID_View edges;
  WhetStone::MFD3D_Diffusion mfd(mesh);
  int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);

  auto p = p_vec.ViewComponent<MirrorHost>("edge", false);
  for (int c = 0; c < ncells; c++) {
    double volume = mesh->cell_volume(c);

    mesh->cell_get_edges(c, edges);
    int nedges = edges.size();
    std::vector<double> cell_solution(nedges);

    for (int k = 0; k < nedges; k++) {
      int e = edges[k];
      cell_solution[k] = p(e,0);

      const AmanziGeometry::Point& xe = mesh->edge_centroid(e);
      double tmp = ana.pressure_exact(xe, t);
      l2_err += std::pow(tmp - p(e,0), 2.0) * volume / nedges;
      inf_err = std::max(inf_err, fabs(tmp - p(e,0)));
      pnorm += std::pow(tmp, 2.0) * volume / nedges;
      // std::cout << e << " at " << xe << " error: " << tmp << " " << p(e,0)
      // << std::endl;
    }
  }
#ifdef HAVE_MPI
  GlobalOp(*mesh->get_comm(), "sum", &pnorm, 1);
  GlobalOp(*mesh->get_comm(), "sum", &l2_err, 1);
  GlobalOp(*mesh->get_comm(), "max", &inf_err, 1);
#endif
  pnorm = sqrt(pnorm);
  l2_err = sqrt(l2_err);
}


/* ******************************************************************
* Error in edge-based fields
****************************************************************** */
inline
void ComputeEdgeMomentsError(const AnalyticBase& ana,
        const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                             const CompositeVector& p_vec, double t, int ngauss, 
                             double& pnorm, double& l2_err, double& inf_err)
{
  pnorm = 0.0;
  l2_err = 0.0;
  inf_err = 0.0;

  AmanziMesh::Entity_ID n0, n1;
  AmanziGeometry::Point x0(ana.dimension()), x1(ana.dimension()), xv(ana.dimension());
  AmanziMesh::Entity_ID_View edges;
  WhetStone::MFD3D_Diffusion mfd(mesh);
  int ncells = mesh->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);

  auto p = p_vec.ViewComponent<MirrorHost>("edge_moments", false); // likely this will fail, not sure what it sould be.  fix when it fails.
  for (int c = 0; c < ncells; c++) {
    double volume = mesh->cell_volume(c);

    mesh->cell_get_edges(c, edges);
    int nedges = edges.size();

    for (int k = 0; k < nedges; k++) {
      int e = edges[k];

      mesh->edge_get_nodes(e, &n0, &n1);
      mesh->node_get_coordinates(n0, &x0);
      mesh->node_get_coordinates(n1, &x1);

      double s0(0.0), s1(0.0);
      for (int n = 0; n <= ngauss; ++n) {
        double gp = WhetStone::q1d_points[ngauss - 1][n];
        double gw = WhetStone::q1d_weights[ngauss - 1][n];

        xv = x0 * gp + x1 * (1.0 - gp);
        s0 += gw * ana.pressure_exact(xv, t);
        s1 += gw * ana.pressure_exact(xv, t) * (0.5 - gp);
      } 

      l2_err += std::pow(s0 - p(e,0), 2.0) * volume / nedges;
      inf_err = std::max(inf_err, fabs(s0 - p(e,0)));
      pnorm += std::pow(s0, 2.0) * volume / nedges;
      // std::cout << e << " at " << (x0 + x1) / 2 << " error: " << s0 << " " <<
      // p(e,0) << std::endl;
    }
  }
#ifdef HAVE_MPI
  GlobalOp(*mesh->get_comm(), "sum", &pnorm, 1);
  GlobalOp(*mesh->get_comm(), "sum", &l2_err, 1);
  GlobalOp(*mesh->get_comm(), "max", &inf_err, 1);
#endif
  pnorm = sqrt(pnorm);
  l2_err = sqrt(l2_err);
}

#endif
