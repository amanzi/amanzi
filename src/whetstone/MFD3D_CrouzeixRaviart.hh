/*
  WhetStone, version 2.1
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Crouzeix-Raviart element: degrees of freedom are moments on faces
  and inside cell.
*/

#ifndef AMANZI_MFD3D_CROUZEIX_RAVIART_HH_
#define AMANZI_MFD3D_CROUZEIX_RAVIART_HH_

#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "Point.hh"

#include "DenseMatrix.hh"
#include "MFD3D.hh"
#include "Tensor.hh"

namespace Amanzi {
namespace WhetStone {

class MFD3D_CrouzeixRaviart : public virtual MFD3D { 
 public:
  MFD3D_CrouzeixRaviart(Teuchos::RCP<const AmanziMesh::Mesh> mesh)
    : MFD3D(mesh),
      InnerProduct(mesh),
      order_(1) {};
  ~MFD3D_CrouzeixRaviart() {};

  // required methods
  // -- mass matrices
  virtual int L2consistency(int c, const Tensor& T, DenseMatrix& N, DenseMatrix& Mc, bool symmetry) { return -1; }
  virtual int MassMatrix(int c, const Tensor& T, DenseMatrix& M) { return -1; } 

  // -- inverse mass matrices
  virtual int L2consistencyInverse(int c, const Tensor& T, DenseMatrix& R, DenseMatrix& Wc, bool symmetry) { return -1; }
  virtual int MassMatrixInverse(int c, const Tensor& T, DenseMatrix& M) { return -1; } 

  // -- stiffness matrix
  virtual int H1consistency(int c, const Tensor& T, DenseMatrix& N, DenseMatrix& Ac);
  virtual int StiffnessMatrix(int c, const Tensor& T, DenseMatrix& A);

  // -- other matrices
  virtual int DivergenceMatrix(int c, DenseMatrix& A) { return -1; }
  virtual int AdvectionMatrix(int c, const std::vector<AmanziGeometry::Point>& u, DenseMatrix& A) { return -1; }

  // -- not relevant or unsupported members
  virtual int MassMatrixPoly(int c, const Polynomial& K, DenseMatrix& M) { return -1; }
  virtual int StiffnessMatrixPoly(int c, const Polynomial& K, DenseMatrix& A) { return -1; }
  virtual int AdvectionMatrix(int c, const AmanziGeometry::Point v, DenseMatrix& A, bool grad_on_test) { return -1; }
  virtual int AdvectionMatrixPoly(int c, const VectorPolynomial& v, DenseMatrix& A, bool grad_on_test) { return -1; }

  // high-order methods
  int H1consistencyHO(int c, int order, const Tensor& T,
                      DenseMatrix& N, DenseMatrix& R, DenseMatrix& Ac, DenseMatrix& G);
  int StiffnessMatrixHO(int c, int order, const Tensor& T, DenseMatrix& A);

  // miscalleneous
  void set_order(int order) { order_ = order; }

 private:
  int order_;
};

}  // namespace WhetStone
}  // namespace Amanzi

#endif

