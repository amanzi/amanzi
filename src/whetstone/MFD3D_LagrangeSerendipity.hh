/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Serendipity Lagrange-type element: degrees of freedom are nodal values 
  and moments on edges, faces and inside cell. The number of later is 
  reduced significantly for polygonal cells. 
*/

#ifndef AMANZI_MFD3D_LAGRANGE_SERENDIPITY_HH_
#define AMANZI_MFD3D_LAGRANGE_SERENDIPITY_HH_

#include "Teuchos_RCP.hpp"

#include "MeshLight.hh"
#include "Point.hh"

#include "BilinearFormFactory.hh"
#include "DenseMatrix.hh"
#include "MFD3D_LagrangeAnyOrder.hh"
#include "Polynomial.hh"
#include "Tensor.hh"

namespace Amanzi {
namespace WhetStone {

class MFD3D_LagrangeSerendipity : public MFD3D_LagrangeAnyOrder { 
 public:
  MFD3D_LagrangeSerendipity(const Teuchos::ParameterList& plist,
                            const Teuchos::RCP<const AmanziMesh::MeshLight>& mesh);

  // required methods
  // -- schema
  virtual std::vector<SchemaItem> schema() const override;

  // -- stiffness matrix
  int H1consistency(int c, const Tensor& T, DenseMatrix& N, DenseMatrix& Ac);
  virtual int StiffnessMatrix(int c, const Tensor& T, DenseMatrix& A) override;

  // -- l2 projectors
  virtual void L2Cell(int c, const std::vector<Polynomial>& ve,
                      const std::vector<Polynomial>& vf,
                      const Polynomial* moments, Polynomial& uc) override {
    ProjectorCell_<AmanziMesh::MeshLight>(mesh_, c, ve, vf, ProjectorType::L2, moments, uc);
  }

  virtual void L2Face(int f, const std::vector<Polynomial>& ve,
                      const Polynomial* moments, Polynomial& uf) override {
    ProjectorFace_(f, ve, ProjectorType::L2, moments, uf);
  }

  // -- h1 projectors
  virtual void H1Cell(int c, const std::vector<Polynomial>& ve,
                      const std::vector<Polynomial>& vf,
                      const Polynomial* moments, Polynomial& uc) override {
    ProjectorCell_<AmanziMesh::MeshLight>(mesh_, c, ve, vf, ProjectorType::H1, moments, uc);
  }

  virtual void H1Face(int f, const std::vector<Polynomial>& ve,
                      const Polynomial* moments, Polynomial& uf) override {
    ProjectorFace_(f, ve, ProjectorType::H1, moments, uf);
  }

  // other methods
  void L2Cell_LeastSquare(int c, const std::vector<Polynomial>& vf,
                          const Polynomial* moments, Polynomial& uc) {
    ProjectorCell_<AmanziMesh::MeshLight>(mesh_, c, vf, vf, ProjectorType::LS, moments, uc);
  }

 private:
  template<class MyMesh>
  void ProjectorCell_(const Teuchos::RCP<const MyMesh>& mymesh, 
                      int c, const std::vector<Polynomial>& ve,
                      const std::vector<Polynomial>& vf,
                      const ProjectorType type,
                      const Polynomial* moments, Polynomial& uc);

  void ProjectorFace_(int f, const std::vector<Polynomial>& ve,
                      const ProjectorType type,
                      const Polynomial* moments, Polynomial& uf);

  template<class MyMesh>
  void CalculateDOFsOnBoundary_(
      const Teuchos::RCP<const MyMesh>& mymesh, 
      int c, const std::vector<Polynomial>& ve,
      const std::vector<Polynomial>& vf, DenseVector& vdof);

 private:
  static RegisteredFactory<MFD3D_LagrangeSerendipity> factory_;
};


/* ******************************************************************
* Generic projector on space of polynomials of order k in cell c.
****************************************************************** */
template<class MyMesh>
void MFD3D_LagrangeSerendipity::ProjectorCell_(
    const Teuchos::RCP<const MyMesh>& mymesh, 
    int c, const std::vector<Polynomial>& ve,
    const std::vector<Polynomial>& vf,
    const ProjectorType type,
    const Polynomial* moments, Polynomial& uc)
{
  // input mesh may have a different dimension than base mesh
  int d = mymesh->space_dimension();

  // selecting regularized basis
  Polynomial ptmp;
  Basis_Regularized<MyMesh> basis;
  basis.Init(mymesh, c, order_, ptmp);

  // calculate stiffness matrix
  Tensor T(d, 1);
  T(0, 0) = 1.0;

  DenseMatrix N, A;
  if (d == 2)
    H1consistency2D_<MyMesh>(mymesh, c, T, N, A);
  else
    H1consistency3D_(c, T, N, A, false);

  // select number of non-aligned edges: we assume cell convexity 
  int nfaces;
  { 
    const auto& faces = mymesh->cell_get_faces(c);
    nfaces = faces.size();
  }
  int eta(3);
  if (nfaces > 3) eta = 4;

  // degrees of freedom: serendipity space S contains all boundary dofs
  // plus a few internal dofs that depend on the value of eta.
  int nd = PolynomialSpaceDimension(d, order_);
  int ndof = N.NumRows();
  int ndof_c = PolynomialSpaceDimension(d, order_ - 2);
  int ndof_cs = PolynomialSpaceDimension(d, order_ - eta);
  int ndof_f(ndof - ndof_c);

  // extract submatrix
  DenseMatrix Ns, NN(nd, nd);
  Ns = N.SubMatrix(0, ndof_f, 0, nd);

  NN.Multiply(Ns, Ns, true);
  NN.InverseSPD();

  // calculate degrees of freedom (Ns^T Ns)^{-1} Ns^T v
  // for consistency with other code, we use v5 for polynomial coefficients
  const AmanziGeometry::Point& xc = mymesh->cell_centroid(c);
  DenseVector v1(nd), v3(std::max(1, ndof_cs)), v5(nd);

  DenseVector vdof(ndof_f + ndof_cs);
  CalculateDOFsOnBoundary_<MyMesh>(mymesh, c, ve, vf, vdof);

  // DOFs inside cell: copy moments from input data
  if (ndof_cs > 0) {
    AMANZI_ASSERT(moments != NULL);
    const DenseVector& v4 = moments->coefs();
    AMANZI_ASSERT(ndof_cs == v4.NumRows());

    for (int n = 0; n < ndof_cs; ++n) {
      vdof(ndof_f + n) = v4(n);
    }
  }

  Ns.Multiply(vdof, v1, true);
  NN.Multiply(v1, v5, false);

  // this gives the least square projector
  uc = basis.CalculatePolynomial(mymesh, c, order_, v5);

  // H1 projector needs to populate moments from ndof_cs + 1 till ndof_c
  if (type == ProjectorType::H1) {
    DenseVector v4(nd);
    DenseMatrix M;
    Polynomial poly(d, order_);

    NumericalIntegration numi(mymesh);
    numi.UpdateMonomialIntegralsCell(c, 2 * order_, integrals_);

    GrammMatrix(poly, integrals_, basis, M);
    M.Multiply(v5, v4, false);

    vdof.Reshape(ndof_f + ndof_c);
    for (int n = ndof_cs; n < ndof_c; ++n) {
      vdof(ndof_f + n) = v4(n) / mymesh->cell_volume(c); 
    }

    R_.Multiply(vdof, v4, true);
    G_.Multiply(v4, v5, false);

    v5(0) = uc(0);
    uc = basis.CalculatePolynomial(mymesh, c, order_, v5);
  }

  // L2 projector is different if the set S contains some internal dofs
  if (type == ProjectorType::L2 && ndof_cs > 0) {
    DenseVector v4(nd), v6(nd - ndof_cs);
    DenseMatrix M, M2;
    Polynomial poly(d, order_);

    NumericalIntegration numi(mymesh);
    numi.UpdateMonomialIntegralsCell(c, 2 * order_, integrals_);

    GrammMatrix(poly, integrals_, basis, M);
    M2 = M.SubMatrix(ndof_cs, nd, 0, nd);
    M2.Multiply(v5, v6, false);

    for (int n = 0; n < ndof_cs; ++n) {
      v4(n) = v3(n) * mymesh->cell_volume(c);
    }

    for (int n = 0; n < nd - ndof_cs; ++n) {
      v4(ndof_cs + n) = v6(n);
    }

    M.InverseSPD();
    M.Multiply(v4, v5, false);

    uc = basis.CalculatePolynomial(mymesh, c, order_, v5);
  }

  // set correct origin 
  uc.set_origin(xc);
}


/* ******************************************************************
* Calculate boundary degrees of freedom in 2D and 3D.
****************************************************************** */
template<class MyMesh>
void MFD3D_LagrangeSerendipity::CalculateDOFsOnBoundary_(
    const Teuchos::RCP<const MyMesh>& mymesh, 
    int c, const std::vector<Polynomial>& ve,
    const std::vector<Polynomial>& vf, DenseVector& vdof)
{
  // input mesh may have a different dimension than base mesh
  int d = mymesh->space_dimension();

  Entity_ID_List nodes, edges;
  mymesh->cell_get_nodes(c, &nodes);
  int nnodes = nodes.size();

  const auto& faces = mymesh->cell_get_faces(c);
  int nfaces = faces.size();

  std::vector<const PolynomialBase*> polys(2);
  NumericalIntegration numi(mymesh);

  int i0, i1, pos;
  AmanziGeometry::Point xv(d);

  // number of moments of faces
  Polynomial pf;
  if (order_ > 1) {
    pf.Reshape(d - 1, order_ - 2);
  }

  int row(nnodes);
  for (int n = 0; n < nfaces; ++n) {
    int f = faces[n];

    Entity_ID_List face_nodes;
    mymesh->face_get_nodes(f, &face_nodes);
    int nfnodes = face_nodes.size();

    if (d == 2) {
      for (int j = 0; j < nfnodes; j++) {
        int v = face_nodes[j];
        mymesh->node_get_coordinates(v, &xv);

        pos = std::distance(nodes.begin(), std::find(nodes.begin(), nodes.end(), v));
        vdof(pos) = vf[n].Value(xv);
      }
    }

    if (order_ > 1) { 
      double area = mymesh->face_area(f);
      const AmanziGeometry::Point& xf = mymesh->face_centroid(f); 
      const AmanziGeometry::Point& normal = mymesh->face_normal(f);

      // local coordinate system with origin at face centroid
      SurfaceCoordinateSystem coordsys(xf, normal);
      const auto& tau = *coordsys.tau();

      polys[0] = &(vf[n]);

      for (auto it = pf.begin(); it < pf.end(); ++it) {
        const int* index = it.multi_index();
        double factor = (d == 2) ? 1.0 : std::pow(area, -(double)it.MonomialSetOrder() / 2);
        Polynomial fmono(d - 1, index, factor);
        fmono.InverseChangeCoordinates(xf, tau);  

        polys[1] = &fmono;

        vdof(row) = numi.IntegratePolynomialsFace(f, polys) / area;
        row++;
      }
    }
  }

  if (d == 3) {
    mymesh->cell_get_edges(c, &edges);
    int nedges = edges.size();

    Polynomial pe(d - 2, order_ - 2);

    for (int n = 0; n < nedges; ++n) {
      int e = edges[n];

      // nodal DOFs
      mymesh->edge_get_nodes(e, &i0, &i1);

      mymesh->node_get_coordinates(i0, &xv);
      pos = std::distance(nodes.begin(), std::find(nodes.begin(), nodes.end(), i0));
      vdof(pos) = ve[n].Value(xv);

      mymesh->node_get_coordinates(i1, &xv);
      pos = std::distance(nodes.begin(), std::find(nodes.begin(), nodes.end(), i1));
      vdof(pos) = ve[n].Value(xv);

      // edge moments
      const auto& xe = mymesh->edge_centroid(e);
      double length = mymesh->edge_length(e);
      std::vector<AmanziGeometry::Point> tau(1, mymesh->edge_vector(e));

      polys[0] = &(ve[n]);

      for (auto it = pe.begin(); it < pe.end(); ++it) {
        const int* index = it.multi_index();
        Polynomial fmono(d - 2, index, 1.0);
        fmono.InverseChangeCoordinates(xe, tau);  

        polys[1] = &fmono;

        vdof(row) = numi.IntegratePolynomialsEdge(e, polys) / length;
        row++;
      }
    }
  }
}

}  // namespace WhetStone
}  // namespace Amanzi

#endif
