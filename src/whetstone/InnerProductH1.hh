/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  The base class for H1 inner products.
*/

#ifndef AMANZI_INNER_PRODUCT_H1_HH_
#define AMANZI_INNER_PRODUCT_H1_HH_

#include "DenseMatrix.hh"
#include "InnerProduct.hh"
#include "Tensor.hh"

namespace Amanzi {
namespace WhetStone {

class Polynomial;

class InnerProductH1 : public virtual InnerProduct { 
 public:
  InnerProductH1() {};
  ~InnerProductH1() {};

  // stiffness matrices
  virtual int H1consistency(int c, const Tensor& T, DenseMatrix& N, DenseMatrix& Mc) = 0;
  virtual int StiffnessMatrix(int c, const Tensor& T, DenseMatrix& A) = 0; 

  // H1 projectors, moments in the optional argument
  virtual void H1Cell(int c, const std::vector<Polynomial>& vf,
                      const Polynomial* moments, Polynomial& vc) {
    Errors::Message msg("H1 projector is not supported/implemented for this scheme.");
    Exceptions::amanzi_throw(msg);
  }
};

}  // namespace WhetStone
}  // namespace Amanzi

#endif

