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

#include "BilinearFormFactory.hh"
#include "DenseMatrix.hh"
#include "DeRham_Edge.hh"
#include "MFD3D.hh"
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

 protected:
  template<typename MyMesh>
  int L2consistency2D_(const Teuchos::RCP<const MyMesh>& mymesh,
                       int c, const Tensor& K, DenseMatrix& N, DenseMatrix& Mc, bool symmetry);

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
    int c, const Tensor& K, DenseMatrix& N, DenseMatrix& Mc, bool symmetry)
{
  // input mesh may have a different dimension than base mesh
  int d = mymesh->space_dimension();

  Entity_ID_List faces;
  mymesh->cell_get_edges(c, &edges);
  int nedges = edges.size();

  int ndc = PolynomialSpaceDimension(d, order_);
  int nde = PolynomialSpaceDimension(d - 1, order_);
  N.Reshape(nedges * nde, ndc * d);

  int row(0);
  for (int i = 0; i < nedges; ++i) {
    int e = edges[i];
    double len = mymesh->edge_area(f);
    const AmanziGeometry::Point& normal = mesh_->edge_vector(e);
      std::vector<AmanziGeometry::Point> tau_edge(1, tau);

      for (auto it = pe.begin(); it < pe.end(); ++it) {
        int m = it.PolynomialPosition();
        const int* index = it.multi_index();
        Polynomial emono(d_, index, 1.0 / len);
        emono.InverseChangeCoordinates(xe, tau_edge);  

        for (auto jt = pc.begin(1); jt < pc.end(); ++jt) {
          const int* index = jt.multi_index();
          double factor = basis.monomial_scales()[jt.MonomialSetOrder()];
          Polynomial cmono(d_, index, factor);
          cmono.set_origin(mesh_->cell_centroid(c));

          polys[0] = &cmono;
          polys[1] = &emono;

          double val = numi.IntegratePolynomialsEdge(e, polys) / len;
          for (int k = 0; k < d_; ++k) Nf(rowf + m, d_ + k) = val * tau[k] / len;
        }
        rowf += nde * d_;
      }
    }

    // lowest-order implementation (for testing only)
    WhetStone::DenseVector v(d_), p0v(nfedges);

    for (int k = 0; k < d_; ++k) {
      AmanziGeometry::Point p(d_);
      p[k] = 1.0;
      auto tmp = ((xf - xc) ^ p) ^ normal;

      for (int l = 0; l < d_; ++l) v(l) = tmp[l];
      P0.Multiply(v, p0v, false);

      for (int i = 0; i < nfedges; ++i) {
        int e = fedges[i];
        int pos = std::distance(edges.begin(), std::find(edges.begin(), edges.end(), e));
        N(pos, k) += p0v(i) * fdirs[n] / 2;
      }
}

}  // namespace WhetStone
}  // namespace Amanzi

#endif
