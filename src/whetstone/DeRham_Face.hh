/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Derham complex: mimetic inner products on faces.
*/

#ifndef AMANZI_DERHAM_FACE_HH_
#define AMANZI_DERHAM_FACE_HH_

#include "Teuchos_RCP.hpp"

#include "MeshLight.hh"

#include "BilinearForm.hh"
#include "DenseMatrix.hh"
#include "MFD3D.hh"
#include "Tensor.hh"

namespace Amanzi {
namespace WhetStone {

class DeRham_Face : public MFD3D {
 public:
  DeRham_Face(const Teuchos::RCP<const AmanziMesh::MeshLight>& mesh) 
    : BilinearForm(mesh) {};

  virtual std::vector<SchemaItem> schema() const override {
    return std::vector<SchemaItem>(1, std::make_tuple(AmanziMesh::FACE, DOF_Type::SCALAR, 1));
  }

  int L2consistency(int c, const Tensor& T, DenseMatrix& N, DenseMatrix& Mc, bool symmetry);
  virtual int MassMatrix(int c, const Tensor& T, DenseMatrix& M) override; 

  int L2consistencyInverse(int c, const Tensor& T, DenseMatrix& R, DenseMatrix& Wc, bool symmetry);
  virtual int MassMatrixInverse(int c, const Tensor& T, DenseMatrix& W) override; 
};

}  // namespace WhetStone
}  // namespace Amanzi

#endif

