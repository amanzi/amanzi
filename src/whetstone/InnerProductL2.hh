/*
  WhetStone, version 2.1
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Base class for L2 mimetic inner products.
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

  // regular polytope
  virtual int L2consistency(int c, const Tensor& T, DenseMatrix& N, DenseMatrix& Mc, bool symmetry) = 0;
  virtual int MassMatrix(int c, const Tensor& T, DenseMatrix& M) = 0; 

  virtual int L2consistencyInverse(int c, const Tensor& T, DenseMatrix& R, DenseMatrix& Wc, bool symmetry) = 0;
  virtual int MassMatrixInverse(int c, const Tensor& T, DenseMatrix& W) = 0; 

  // generalized polytope
  virtual int L2consistencyGeneralized(int c, const Tensor& K, DenseMatrix& N, DenseMatrix& Mc, bool symmetry) = 0;
  virtual int MassMatrixGeneralized(int c, const Tensor& K, DenseMatrix& M) = 0;

  virtual int L2consistencyInverseGeneralized(int c, const Tensor& K, DenseMatrix& N, DenseMatrix& Mc, bool symmetry) = 0;
  virtual int MassMatrixInverseGeneralized(int c, const Tensor& K, DenseMatrix& M) = 0;
};

}  // namespace WhetStone
}  // namespace Amanzi

#endif

