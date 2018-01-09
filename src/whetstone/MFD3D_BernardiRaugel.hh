/*
  WhetStone, version 2.1
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Bernardi Raugel element: velocity space vectors at nodes and 
  normal components on faces.
*/

#ifndef AMANZI_MFD3D_BERNARDI_RAUGEL_HH_
#define AMANZI_MFD3D_BERNARDI_RAUGEL_HH_

#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "Point.hh"

#include "DenseMatrix.hh"
#include "MFD3D.hh"
#include "Tensor.hh"

namespace Amanzi {
namespace WhetStone {

class MFD3D_BernardiRaugel : public virtual MFD3D { 
 public:
  MFD3D_BernardiRaugel(Teuchos::RCP<const AmanziMesh::Mesh> mesh)
    : MFD3D(mesh),
      InnerProduct(mesh) {};
  ~MFD3D_BernardiRaugel() {};

  // required methods
  // -- mass matrices
  virtual int L2consistency(int c, const Tensor& T, DenseMatrix& N, DenseMatrix& Mc, bool symmetry) { return -1; }
  virtual int MassMatrix(int c, const Tensor& T, DenseMatrix& M) { return -1; } 

  // -- inverse mass matrices
  virtual int L2consistencyInverse(int c, const Tensor& T, DenseMatrix& R, DenseMatrix& Wc, bool symmetry) { return -1; }
  virtual int MassMatrixInverse(int c, const Tensor& T, DenseMatrix& M) { return -1; } 

  // -- stiffness matrix
  virtual int H1consistency(int c, const Tensor& T, DenseMatrix& N, DenseMatrix& Mc);
  virtual int StiffnessMatrix(int c, const Tensor& T, DenseMatrix& A);

  // -- other matrices
  virtual int DivergenceMatrix(int c, DenseMatrix& A);
  virtual int AdvectionMatrix(int c, const std::vector<AmanziGeometry::Point>& u, DenseMatrix& A);

  // -- not relevant or unsupported members
  virtual int MassMatrixPoly(int c, const Polynomial& K, DenseMatrix& M) { return -1; }
  virtual int StiffnessMatrixPoly(int c, const Polynomial& K, DenseMatrix& A) { return -1; }
  virtual int AdvectionMatrix(int c, const AmanziGeometry::Point v, DenseMatrix& A, bool grad_on_test) { return -1; }
  virtual int AdvectionMatrixPoly(int c, const VectorPolynomial& v, DenseMatrix& A, bool grad_on_test) { return -1; }
};

}  // namespace WhetStone
}  // namespace Amanzi

#endif

