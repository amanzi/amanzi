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
#include "DenseMatrix.hh"
#include "DeRham_Edge.hh"
#include "MFD3D.hh"
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
  virtual int StiffnessMatrix(int c, const Tensor& T, DenseMatrix& A) override { return 0; }

  // access
  PolynomialOnMesh& integrals() { return integrals_; }

 protected:
  template<typename MyMesh>
  int L2consistency2D_(const Teuchos::RCP<const MyMesh>& mymesh,
                       int c, const Tensor& K, DenseMatrix& N, DenseMatrix& Mc, DenseMatrix& MG);

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

  std::vector<const PolynomialBase*> polys(2);

  // iterators
  VectorPolynomialIterator it0(d, d, order_), it1(d, d, order_);
  it0.begin();
  it1.end();

  for (auto it = it0; it < it1; ++it) {
    int k = it.VectorComponent();
    int n = it.VectorPolynomialPosition();
    int m = it.PolynomialPosition();

    const int* index = it.multi_index();
    double factor = basis.monomial_scales()[it.MonomialSetOrder()];
    Polynomial cmono(d, index, factor);
    cmono.set_origin(mymesh->cell_centroid(c));

    int row(0);
    for (int i = 0; i < nedges; ++i) {
      int e = edges[i];
      double len = mymesh->edge_length(e);
      const AmanziGeometry::Point& xe = mymesh->edge_centroid(e);

      const AmanziGeometry::Point& tau = mymesh->edge_vector(e);
      std::vector<AmanziGeometry::Point> tau_edge(1, tau);

      for (auto jt = pe.begin(); jt < pe.end(); ++jt) {
        const int* jndex = jt.multi_index();
        Polynomial emono(d, jndex, tau[k] / len);
        emono.InverseChangeCoordinates(xe, tau_edge);  

        polys[0] = &cmono;
        polys[1] = &emono;

        N(row, n) = numi.IntegratePolynomialsEdge(e, polys) / len;
      }
      row++;
    }
  }

  // calculate Mc = P0 M_G P0^T
  GrammMatrix(numi, order_, integrals_, basis, MG);
  
  DenseMatrix NT, P0T;
  Tensor Id(d, 2);
  Id.MakeDiagonal(1.0);

  NT.Transpose(N);
  auto NN = NT * N;
  NN.InverseMoorePenrose();

  Tensor Kinv(K);
  Kinv.Inverse();

  auto P0 = N * NN;
  P0T.Transpose(P0);
  Mc = P0 * ((Id * Kinv) ^ MG) * P0T;

  return WHETSTONE_ELEMENTAL_MATRIX_OK;
}

}  // namespace WhetStone
}  // namespace Amanzi

#endif
