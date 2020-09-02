/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Derham complex: mimetic inner products on edges.
*/

#ifndef AMANZI_DERHAM_EDGE_HH_
#define AMANZI_DERHAM_EDGE_HH_

#include "Teuchos_RCP.hpp"

#include "Mesh.hh"

#include "BilinearForm.hh"
#include "DenseMatrix.hh"
#include "MFD3D.hh"
#include "Tensor.hh"

namespace Amanzi {
namespace WhetStone {

class DeRham_Edge : public MFD3D {
 public:
  DeRham_Edge(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) 
    : BilinearForm(mesh) {};

  virtual std::vector<SchemaItem> schema() const override {
    return std::vector<SchemaItem>(1, std::make_tuple(AmanziMesh::EDGE, DOF_Type::SCALAR, 1));
  }

  int L2consistency(int c, const Tensor& T, DenseMatrix& N, DenseMatrix& Mc, bool symmetry);
  virtual int MassMatrix(int c, const Tensor& T, DenseMatrix& M) override; 

  int L2consistencyInverse(int c, const Tensor& T, DenseMatrix& R, DenseMatrix& Wc, bool symmetry);
  virtual int MassMatrixInverse(int c, const Tensor& T, DenseMatrix& W) override; 

 protected:
  int L2consistency2D_(int c, const Tensor& T, DenseMatrix& N, DenseMatrix& Mc);
  int L2consistency3D_(int c, const Tensor& T, DenseMatrix& N, DenseMatrix& Mc);

  int L2consistencyInverse2D_(int c, const Tensor& T, DenseMatrix& R, DenseMatrix& Wc);
  int L2consistencyInverse3D_(int c, const Tensor& T, DenseMatrix& R, DenseMatrix& Wc);
};

}  // namespace WhetStone
}  // namespace Amanzi

#endif

