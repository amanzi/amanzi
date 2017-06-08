/*
  WhetStone, version 2.1
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  The base class for mimetic bilinear forms.
*/

#ifndef AMANZI_INNER_PRODUCT_H1_HH_
#define AMANZI_INNER_PRODUCT_H1_HH_

#include "Teuchos_RCP.hpp"

#include "Mesh.hh"

#include "DenseMatrix.hh"
#include "InnerProduct.hh"
#include "Tensor.hh"

namespace Amanzi {
namespace WhetStone {

class InnerProductH1 : virtual InnerProduct { 
 public:
  explicit InnerProductH1(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) : InnerProduct(mesh) {};
  ~InnerProductH1() {};

  virtual int H1consistency(int c, const Tensor& T, DenseMatrix& N, DenseMatrix& Mc) = 0;
  virtual int StiffnessMatrix(int c, const Tensor& T, DenseMatrix& A) = 0; 
};

}  // namespace WhetStone
}  // namespace Amanzi

#endif

