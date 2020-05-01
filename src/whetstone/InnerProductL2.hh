/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  The base class for L2 inner products.
*/

#ifndef AMANZI_INNER_PRODUCT_L2_HH_
#define AMANZI_INNER_PRODUCT_L2_HH_

#include "DenseMatrix.hh"
#include "InnerProduct.hh"
#include "Tensor.hh"

namespace Amanzi {
namespace WhetStone {

class InnerProductL2 : public virtual InnerProduct { 
 public:
  InnerProductL2() {};
  ~InnerProductL2() {};

  // mass matrices
  virtual int L2consistency(int c, const Tensor<>& T, DenseMatrix<>& N, DenseMatrix<>& Mc, bool symmetry) = 0;
  virtual int MassMatrix(int c, const Tensor<>& T, DenseMatrix<>& M) = 0;

  // optional inverse mass matrices
  virtual int L2consistencyInverse(int c, const Tensor<>& T, DenseMatrix<>& R, DenseMatrix<>& Wc, bool symmetry) {
    Errors::Message msg("L2 inverse consistency is not supported for this space.");
    Exceptions::amanzi_throw(msg);
    return WHETSTONE_ELEMENTAL_MATRIX_OK;
  }
  virtual int MassMatrixInverse(int c, const Tensor<>& T, DenseMatrix<>& W) {
    MassMatrix(c, T, W);
    W.Inverse();
    return WHETSTONE_ELEMENTAL_MATRIX_OK;
  }
};

}  // namespace WhetStone
}  // namespace Amanzi

#endif
