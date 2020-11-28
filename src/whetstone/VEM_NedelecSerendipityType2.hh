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
#include "MFD3D.hh"
#include "Monomial.hh"
#include "NumericalIntegration.hh"
#include "Polynomial.hh"
#include "PolynomialOnMesh.hh"
#include "Tensor.hh"

namespace Amanzi {
namespace WhetStone {

class VEM_NedelecSerendipityType2 : public DeRham_Edge { 
 public:
  VEM_NedelecSerendipityType2(const Teuchos::ParameterList& plist,
                              const Teuchos::RCP<const AmanziMesh::MeshLight>& mesh);
  ~VEM_NedelecSerendipityType2() {};

  // required methods
  // -- schema
  virtual std::vector<SchemaItem> schema() const override;

  // -- mass matrix
  int L2consistency(int c, const Tensor& K, DenseMatrix& N, DenseMatrix& Mc, bool symmetry);
  virtual int MassMatrix(int c, const Tensor& K, DenseMatrix& A) override;

  int MassMatrixFace(int f, const Tensor& K, DenseMatrix& M);

  // -- stiffness matrix
  virtual int StiffnessMatrix(int c, const Tensor& T, DenseMatrix& A) override;

  // other methods
  int StiffnessMatrix(int c, const Tensor& K, DenseMatrix& A, DenseMatrix& M, DenseMatrix& C);
  void CurlMatrix(int c, DenseMatrix& C);

  // -- l2 projector
  void L2Cell(int c, const std::vector<VectorPolynomial>& ve,
              const std::vector<VectorPolynomial>& vf,
              const Polynomial* moments, VectorPolynomial& uc) {
    ProjectorCell_(mesh_, c, ve, vf, ProjectorType::L2, moments, uc);
  }

  void L2Face(int f, const std::vector<VectorPolynomial>& ve,
              const Polynomial* moments, VectorPolynomial& uf) {
    ProjectorFace_(f, ve, ProjectorType::L2, moments, uf);
  }

  // access
  PolynomialOnMesh& integrals() { return integrals_; }

  // support
  void CalculateDOFsOnBoundary(
      const Teuchos::RCP<const AmanziMesh::MeshLight>& mymesh, 
      int c, const std::vector<VectorPolynomial>& ve,
      const std::vector<VectorPolynomial>& vf, DenseVector& vdof);

 protected:
  int L2consistency2D_(const Teuchos::RCP<const AmanziMesh::MeshLight>& mymesh,
                       int c, const Tensor& K, DenseMatrix& N, DenseMatrix& Mc, DenseMatrix& MG);

 private:
  void ProjectorCell_(const Teuchos::RCP<const AmanziMesh::MeshLight>& mymesh, 
                      int c, const std::vector<VectorPolynomial>& ve,
                      const std::vector<VectorPolynomial>& vf,
                      const ProjectorType type,
                      const Polynomial* moments, VectorPolynomial& uc);

  void ProjectorFace_(int f, const std::vector<VectorPolynomial>& ve,
                      const ProjectorType type,
                      const Polynomial* moments, VectorPolynomial& uf);

  // auxiliary functions
  void MatrixOfDofs_(
      int c, const Entity_ID_List& edges,
      const Basis_Regularized& basis,
      const NumericalIntegration& numi,
      DenseMatrix& N);

  void L2ProjectorsOnFaces_(
      int c, const Tensor& K, const Entity_ID_List& faces,
      std::vector<WhetStone::DenseMatrix>& vL2f, 
      std::vector<WhetStone::DenseMatrix>& vMGf,
      std::vector<Basis_Regularized>& vbasisf,
      std::vector<std::shared_ptr<WhetStone::SurfaceCoordinateSystem> >& vcoordsys,
      int MGorder);

  void L2ProjectorOnEdge_(WhetStone::DenseMatrix& L2e, int order);

  void L2consistency3DFace_Method1_(
      const VectorPolynomial& p1,
      const VectorPolynomial& xyz,
      const SurfaceCoordinateSystem& coordsys,
      const Basis_Regularized& basis,
      const DenseMatrix& L2f,
      const DenseMatrix& MGf,
      DenseVector& p0v);

  void L2consistency3DFace_Method2_(
      int f,
      const VectorPolynomial& p1,
      const SurfaceCoordinateSystem& coordsys,
      const Basis_Regularized& basis,
      const DenseMatrix& L2f,
      const DenseMatrix& MGf,
      DenseVector& p0v);

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
inline
int VEM_NedelecSerendipityType2::L2consistency2D_(
    const Teuchos::RCP<const AmanziMesh::MeshLight>& mymesh,
    int c, const Tensor& K, DenseMatrix& N, DenseMatrix& Mc, DenseMatrix& MG)
{
  // input mesh may have a different dimension than base mesh
  int d = mymesh->space_dimension();

  const auto& edges = mymesh->cell_get_edges(c);
  int nedges = edges.size();

  Polynomial poly(d, order_);

  int ndc = PolynomialSpaceDimension(d, order_);
  int nde = PolynomialSpaceDimension(d - 1, order_);
  N.Reshape(nedges * nde, ndc * d);

  // selecting regularized basis (parameter integrals is not used)
  PolynomialOnMesh integrals_f;
  integrals_f.set_id(c);

  Basis_Regularized basis;
  basis.Init(mymesh, c, order_, integrals_f.poly());

  // pre-calculate integrals of monomials 
  NumericalIntegration numi(mymesh);
  numi.UpdateMonomialIntegralsCell(c, 2 * order_, integrals_f);

  // iterators
  std::vector<double > moments;
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

      numi.CalculatePolynomialMomentsEdge(e, cmono, order_, moments);
      for (int l = 0; l < moments.size(); ++l) {
        N(row, n) = moments[l] * tau[k] / len;
        row++;
      }
    }
  }

  // calculate Mc = P0 M_G P0^T
  GrammMatrix(numi, order_, integrals_f, basis, MG);
  
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

  return 0;
}


/* ******************************************************************
* Generic projector on space of polynomials of order k in cell c.
****************************************************************** */
inline
void VEM_NedelecSerendipityType2::ProjectorCell_(
    const Teuchos::RCP<const AmanziMesh::MeshLight>& mymesh, 
    int c, const std::vector<VectorPolynomial>& ve,
    const std::vector<VectorPolynomial>& vf,
    const ProjectorType type,
    const Polynomial* moments, VectorPolynomial& uc)
{
  // input mesh may have a different dimension than base mesh
  int d = mymesh->space_dimension();

  // selecting regularized basis
  Polynomial ptmp;
  Basis_Regularized basis;
  basis.Init(mymesh, c, order_, ptmp);

  // calculate stiffness matrix
  Tensor T(d, 1);
  T(0, 0) = 1.0;

  DenseMatrix N, Mc, MG;
  if (d == 2) 
    L2consistency2D_(mymesh, c, T, N, Mc, MG);
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
  CalculateDOFsOnBoundary(mymesh, c, ve, vf, vdof);

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
inline
void VEM_NedelecSerendipityType2::CalculateDOFsOnBoundary(
    const Teuchos::RCP<const AmanziMesh::MeshLight>& mymesh, 
    int c, const std::vector<VectorPolynomial>& ve,
    const std::vector<VectorPolynomial>& vf, DenseVector& vdof)
{
  const auto& edges = mymesh->cell_get_edges(c);
  int nedges = edges.size();

  std::vector<const WhetStoneFunction*> funcs(2);
  NumericalIntegration numi(mymesh);

  // number of moments on edges
  std::vector<double> moments;

  int row(0);
  for (int n = 0; n < nedges; ++n) {
    int e = edges[n];
    double length = mymesh->edge_length(e);
    const auto& tau = mymesh->edge_vector(e);

    auto poly = ve[n] * tau;

    numi.CalculatePolynomialMomentsEdge(e, poly, order_, moments);
    for (int k = 0; k < moments.size(); ++k) {
      vdof(row) = moments[k] / length;
      row++;
    }
  }
}

}  // namespace WhetStone
}  // namespace Amanzi

#endif
