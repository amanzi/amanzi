/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  High-order 3D Nedelec serendipity element of type 2: degrees of 
  freedom are moments on edges, faces and incide cell.
*/

#ifndef AMANZI_VEM_NEDELEC_SERENDIPITY_TYPE2_HH_
#define AMANZI_VEM_NEDELEC_SERENDIPITY_TYPE2_HH_

#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "Point.hh"

#include "Basis_Regularized.hh"
#include "BilinearFormFactory.hh"
#include "GrammMatrix.hh"
#include "DenseMatrix.hh"
#include "DeRham_Edge.hh"
#include "FunctionPower.hh"
#include "MFD3D.hh"
#include "Monomial.hh"
#include "NumericalIntegration.hh"
#include "Polynomial.hh"
#include "PolynomialOnMesh.hh"
#include "SurfaceMiniMesh.hh"
#include "Tensor.hh"

namespace Amanzi {
namespace WhetStone {

class VEM_NedelecSerendipityType2 : public MFD3D,
                                    public DeRham_Edge { 
 public:
  VEM_NedelecSerendipityType2(const Teuchos::ParameterList& plist,
                            const Teuchos::RCP<const AmanziMesh::Mesh>& mesh);
  ~VEM_NedelecSerendipityType2() {};

  // required methods
  // -- schema
  virtual std::vector<SchemaItem> schema() const override;

  // -- mass matrix
  virtual int L2consistency(int c, const Tensor& K, DenseMatrix& N, DenseMatrix& Mc, bool symmetry) override;
  virtual int MassMatrix(int c, const Tensor& K, DenseMatrix& A) override;

  // -- stiffness matrix
  virtual int H1consistency(int c, const Tensor& T, DenseMatrix& N, DenseMatrix& Ac) override { return 0; }
  virtual int StiffnessMatrix(int c, const Tensor& T, DenseMatrix& A) override;

  // other methods
  void CurlMatrix(int c, DenseMatrix& C);

  // -- l2 projector
  void L2Cell(int c, const std::vector<VectorPolynomial>& ve,
              const std::vector<VectorPolynomial>& vf,
              const Polynomial* moments, VectorPolynomial& uc) {
    ProjectorCell_<AmanziMesh::Mesh>(mesh_, c, ve, vf, ProjectorType::L2, moments, uc);
  }

  void L2Face(int f, const std::vector<VectorPolynomial>& ve,
              const Polynomial* moments, VectorPolynomial& uf) {
    ProjectorFace_(f, ve, ProjectorType::L2, moments, uf);
  }

  // access
  PolynomialOnMesh& integrals() { return integrals_; }

  // support
  template<class MyMesh>
  void CalculateDOFsOnBoundary(
      const Teuchos::RCP<const MyMesh>& mymesh, 
      int c, const std::vector<VectorPolynomial>& ve,
      const std::vector<VectorPolynomial>& vf, DenseVector& vdof);

 protected:
  template<typename MyMesh>
  int L2consistency2D_(const Teuchos::RCP<const MyMesh>& mymesh,
                       int c, const Tensor& K, DenseMatrix& N, DenseMatrix& Mc, DenseMatrix& MG);

 private:
  template<class MyMesh>
  void ProjectorCell_(const Teuchos::RCP<const MyMesh>& mymesh, 
                      int c, const std::vector<VectorPolynomial>& ve,
                      const std::vector<VectorPolynomial>& vf,
                      const ProjectorType type,
                      const Polynomial* moments, VectorPolynomial& uc);

  void ProjectorFace_(int f, const std::vector<VectorPolynomial>& ve,
                      const ProjectorType type,
                      const Polynomial* moments, VectorPolynomial& uf);

 protected:
  using MFD3D::mesh_;
  using MFD3D::d_;

 private:
  PolynomialOnMesh integrals_;

 private:
  static RegisteredFactory<VEM_NedelecSerendipityType2> factory_;
};


/* ******************************************************************
* High-order consistency condition for the 2D mass matrix. 
****************************************************************** */
template <class MyMesh>
int VEM_NedelecSerendipityType2::L2consistency2D_(
    const Teuchos::RCP<const MyMesh>& mymesh,
    int c, const Tensor& K, DenseMatrix& N, DenseMatrix& Mc, DenseMatrix& MG)
{
  // input mesh may have a different dimension than base mesh
  int d = mymesh->space_dimension();

  Entity_ID_List edges;
  mymesh->cell_get_edges(c, &edges);
  int nedges = edges.size();

  Polynomial poly(d, order_), pe(d - 1, order_);

  int ndc = PolynomialSpaceDimension(d, order_);
  int nde = PolynomialSpaceDimension(d - 1, order_);
  N.Reshape(nedges * nde, ndc * d);

  // selecting regularized basis (parameter integrals is not used)
  Basis_Regularized<MyMesh> basis;
  basis.Init(mymesh, c, order_, integrals_.poly());

  // pre-calculate integrals of monomials 
  NumericalIntegration<MyMesh> numi(mymesh);
  numi.UpdateMonomialIntegralsCell(c, 2 * order_, integrals_);

  std::vector<const WhetStoneFunction*> funcs(2);

  // iterators
  VectorPolynomialIterator it0(d, d, order_), it1(d, d, order_);
  it0.begin();
  it1.end();

  for (auto it = it0; it < it1; ++it) {
    int k = it.VectorComponent();
    int n = it.VectorPolynomialPosition();
    int m = it.MonomialSetOrder();

    const int* index = it.multi_index();
    double factor = basis.monomial_scales()[m];
    Monomial cmono(d, index, factor);
    cmono.set_origin(mymesh->cell_centroid(c));

    int row(0);
    for (int i = 0; i < nedges; ++i) {
      int e = edges[i];
      double len = mymesh->edge_length(e);
      const AmanziGeometry::Point& tau = mymesh->edge_vector(e);
      std::vector<AmanziGeometry::Point> tau_edge(1, tau);

      for (auto jt = pe.begin(); jt < pe.end(); ++jt) {
        int m2 = jt.MonomialSetOrder();
        FunctionPower efunc(tau[k] / len, m2);

        funcs[0] = &cmono;
        funcs[1] = &efunc;

        N(row, n) = numi.IntegrateFunctionsEdge(e, funcs, m + m2) / len;
        row++;
      }
    }
  }

  // calculate Mc = P0 M_G P0^T
  GrammMatrix(numi, order_, integrals_, basis, MG);
  
  DenseMatrix NT, P0T;
  Tensor Id(d, 2);
  Id.MakeDiagonal(1.0);

  NT.Transpose(N);
  auto NN = NT * N;
  NN.InverseSPD();

  Tensor Kinv(K);
  Kinv.Inverse();

  auto P0 = N * NN;
  P0T.Transpose(P0);
  Mc = P0 * ((Id * Kinv) ^ MG) * P0T;

  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}


/* ******************************************************************
* Generic projector on space of polynomials of order k in cell c.
****************************************************************** */
template<class MyMesh>
void VEM_NedelecSerendipityType2::ProjectorCell_(
    const Teuchos::RCP<const MyMesh>& mymesh, 
    int c, const std::vector<VectorPolynomial>& ve,
    const std::vector<VectorPolynomial>& vf,
    const ProjectorType type,
    const Polynomial* moments, VectorPolynomial& uc)
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

  DenseMatrix N, Mc, MG;
  if (d == 2) 
    L2consistency2D_<MyMesh>(mymesh, c, T, N, Mc, MG);
  else
    L2consistency(c, T, N, Mc, true);

  // select number of non-aligned faces: we assume cell convexity 
  // int eta(2);

  // degrees of freedom: serendipity space S contains all boundary dofs
  // plus a few internal dofs that depend on the value of eta.
  int ndof = N.NumRows();
  int ndof_cs = 0;  // required cell moments
  int ndof_s(ndof);  // serendipity dofs

  // extract submatrix
  int ncols = N.NumCols();
  DenseMatrix Ns, NN(ncols, ncols);
  Ns = N.SubMatrix(0, ndof_s, 0, ncols);

  NN.Multiply(Ns, Ns, true);
  NN.InverseSPD();

  // calculate degrees of freedom (Ns^T Ns)^{-1} Ns^T v
  // for consistency with other code, we use v5 for polynomial coefficients
  const AmanziGeometry::Point& xc = mymesh->cell_centroid(c);
  DenseVector v1(ncols), v5(ncols);

  DenseVector vdof(ndof_s + ndof_cs);
  CalculateDOFsOnBoundary<MyMesh>(mymesh, c, ve, vf, vdof);

  // DOFs inside cell: copy moments from input data
  if (ndof_cs > 0) {
    AMANZI_ASSERT(moments != NULL);
    const DenseVector& v3 = moments->coefs();
    AMANZI_ASSERT(ndof_cs == v3.NumRows());

    for (int n = 0; n < ndof_cs; ++n) {
      vdof(ndof_s + n) = v3(n);
    }
  }

  Ns.Multiply(vdof, v1, true);
  NN.Multiply(v1, v5, false);

  // this gives the least square projector
  int stride = v5.NumRows() / d;
  DenseVector v4(stride);

  uc.resize(d);
  for (int k = 0; k < d; ++k) {
    for (int i = 0; i < stride; ++i) v4(i) = v5(k * stride + i);
    uc[k] = basis.CalculatePolynomial(mymesh, c, order_, v4);
  }

  // set correct origin 
  uc.set_origin(xc);
}


/* ******************************************************************
* Calculate boundary degrees of freedom in 2D and 3D.
****************************************************************** */
template<class MyMesh>
void VEM_NedelecSerendipityType2::CalculateDOFsOnBoundary(
    const Teuchos::RCP<const MyMesh>& mymesh, 
    int c, const std::vector<VectorPolynomial>& ve,
    const std::vector<VectorPolynomial>& vf, DenseVector& vdof)
{
  Entity_ID_List edges;
  mymesh->cell_get_edges(c, &edges);
  int nedges = edges.size();

  std::vector<const WhetStoneFunction*> funcs(2);
  NumericalIntegration<MyMesh> numi(mymesh);

  // number of moments on edges
  Polynomial pe(1, order_);

  int row(0);
  for (int n = 0; n < nedges; ++n) {
    int e = edges[n];
    double length = mymesh->edge_length(e);
    const auto& tau = mymesh->edge_vector(e);

    auto poly = ve[n] * tau;
    int m = poly.order();

    for (auto it = pe.begin(); it < pe.end(); ++it) {
      int m2 = it.MonomialSetOrder();
      FunctionPower efunc(1.0 / length, m2);

      funcs[0] = &poly;
      funcs[1] = &efunc;

      vdof(row) = numi.IntegrateFunctionsEdge(e, funcs, m + m2) / length;
      row++;
    }
  }
}

}  // namespace WhetStone
}  // namespace Amanzi

#endif
