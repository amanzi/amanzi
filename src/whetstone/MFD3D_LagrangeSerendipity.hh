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

#include "Mesh.hh"
#include "Point.hh"

#include "BilinearFormFactory.hh"
#include "DenseMatrix.hh"
#include "MFD3D_LagrangeAnyOrder.hh"
#include "Polynomial.hh"
#include "SurfaceMiniMesh.hh"
#include "Tensor.hh"

namespace Amanzi {
namespace WhetStone {

class MFD3D_LagrangeSerendipity : public MFD3D_LagrangeAnyOrder { 
 public:
  MFD3D_LagrangeSerendipity(const Teuchos::ParameterList& plist,
                            const Teuchos::RCP<const AmanziMesh::Mesh>& mesh);
  ~MFD3D_LagrangeSerendipity() {};

  // required methods
  // -- schema
  virtual std::vector<SchemaItem> schema() const override;

  // -- stiffness matrix
  virtual int H1consistency(int c, const Tensor& T, DenseMatrix& N, DenseMatrix& Ac) override;
  virtual int StiffnessMatrix(int c, const Tensor& T, DenseMatrix& A) override;

  // -- l2 projectors
  virtual void L2Cell(int c, const std::vector<Polynomial>& ve,
                      const std::vector<Polynomial>& vf,
                      const Polynomial* moments, Polynomial& uc) override {
    ProjectorCell_<AmanziMesh::Mesh>(mesh_, c, ve, vf, ProjectorType::L2, moments, uc);
  }

  // -- h1 projectors
  virtual void H1Cell(int c, const std::vector<Polynomial>& ve,
                      const std::vector<Polynomial>& vf,
                      const Polynomial* moments, Polynomial& uc) override {
    ProjectorCell_<AmanziMesh::Mesh>(mesh_, c, ve, vf, ProjectorType::H1, moments, uc);
  }

  virtual void H1Face(int f, const std::vector<Polynomial>& ve,
                      const Polynomial* moments, Polynomial& uf) override {
    auto coordsys = std::make_shared<SurfaceCoordinateSystem>(mesh_->face_normal(f));
    Teuchos::RCP<const SurfaceMiniMesh> surf_mesh = Teuchos::rcp(new SurfaceMiniMesh(mesh_, coordsys));
    ProjectorCell_<SurfaceMiniMesh>(surf_mesh, f, ve, ve, ProjectorType::H1, moments, uf);
  }

  // other methods
  void L2Cell_LeastSquare(int c, const std::vector<Polynomial>& vf,
                          const Polynomial* moments, Polynomial& uc) {
    ProjectorCell_<AmanziMesh::Mesh>(mesh_, c, vf, vf, ProjectorType::LS, moments, uc);
  }

 private:
  template<class MyMesh>
  void ProjectorCell_(const Teuchos::RCP<const MyMesh>& mymesh, 
                      int c, const std::vector<Polynomial>& ve,
                      const std::vector<Polynomial>& vf,
                      const ProjectorType type,
                      const Polynomial* moments, Polynomial& uc);

  void CalculateDOFsOnBoundary_(
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
  // selecting regularized basis
  Polynomial ptmp;
  Basis_Regularized<MyMesh> basis;
  basis.Init(mymesh, c, order_, ptmp);

  // calculate stiffness matrix
  Tensor T(d_, 1);
  DenseMatrix N, A;

  T(0, 0) = 1.0;
  MFD3D_LagrangeAnyOrder::H1consistency(c, T, N, A);  

  // select number of non-aligned edges: we assume cell convexity 
  int nfaces;
  { 
    Entity_ID_List faces;
    mymesh->cell_get_faces(c, &faces);
    nfaces = faces.size();
  }
  int eta(3);
  if (nfaces > 3) eta = 4;

  // degrees of freedom: serendipity space S contains all boundary dofs
  // plus a few internal dofs that depend on the value of eta.
  int nd = PolynomialSpaceDimension(d_, order_);
  int ndof = A.NumRows();
  int ndof_c = PolynomialSpaceDimension(d_, order_ - 2);
  int ndof_cs = PolynomialSpaceDimension(d_, order_ - eta);
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
  CalculateDOFsOnBoundary_(c, ve, vf, vdof);

  // DOFs inside cell: copy moments from input data
  if (ndof_cs > 0) {
    AMANZI_ASSERT(moments != NULL);
    const DenseVector& v3 = moments->coefs();
    AMANZI_ASSERT(ndof_cs == v3.NumRows());

    for (int n = 0; n < ndof_cs; ++n) {
      vdof(ndof_f + n) = v3(n);
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
    Polynomial poly(d_, order_);

    NumericalIntegration<MyMesh> numi(mymesh);
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
    Polynomial poly(d_, order_);

    NumericalIntegration<MyMesh> numi(mymesh);
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

}  // namespace WhetStone
}  // namespace Amanzi

#endif
