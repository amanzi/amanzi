/*
  WhetStone, version 2.1
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  The base virtual class for factory of mimetic, DG and other 
  schemes on polytopal meshes.
*/

#ifndef AMANZI_BILINEAR_FORM_HH_
#define AMANZI_BILINEAR_FORM_HH_

#include "Teuchos_RCP.hpp"

#include "errors.hh"
#include "Mesh.hh"
#include "Point.hh"

#include "DenseMatrix.hh"
#include "InnerProductL2.hh"
#include "InnerProductH1.hh"
#include "Projectors.hh"
#include "Tensor.hh"
#include "VectorPolynomial.hh"

namespace Amanzi {
namespace WhetStone {

class Polynomial;

class BilinearForm : public virtual InnerProductL2,
                     public virtual InnerProductH1,
                     public Projectors {
 public:
  explicit BilinearForm() : order_(1) {};
  ~BilinearForm() {};

  // additional members
  // -- low-order schemes require constant vector/tensor coefficients
  //    we also specify function to which gradient operator is applied
  virtual int AdvectionMatrix(int c, const AmanziGeometry::Point v, DenseMatrix& A, bool grad_on_test) {;
    Errors::Message msg("AdvectionMatrix: scalar velocity is not supported.");
    Exceptions::amanzi_throw(msg);
    return 1;
  }

  virtual int DivergenceMatrix(int c, DenseMatrix& A) {
    Errors::Message msg("DivergenceMatrix is not supported.");
    Exceptions::amanzi_throw(msg);
    return 1;
  }

  // -- high-order schemes may require polynomial coefficients
  //    We use ending "Poly" to simplify logic of derived classes (no keyword "using")
  virtual int MassMatrixPoly(int c, const Polynomial& K, DenseMatrix& M) {
    Errors::Message msg("MassMatrix: polynomial coefficient is not supported.");
    Exceptions::amanzi_throw(msg);
    return 1;
  }

  virtual int StiffnessMatrixPoly(int c, const Polynomial& K, DenseMatrix& A) {
    Errors::Message msg("StiffnessMatrix: polynomial coefficient is not supported.");
    Exceptions::amanzi_throw(msg);
    return 1;
  }

  virtual int AdvectionMatrixPoly(int c, const VectorPolynomial& v, DenseMatrix& A, bool grad_on_test) {
    Errors::Message msg("AdvectionMatrix: polynomial coefficient is not supported.");
    Exceptions::amanzi_throw(msg);
    return 1;
  }

  // miscalleneous
  void set_order(int order) { order_ = order; }

 protected:
  int order_;
};

}  // namespace WhetStone
}  // namespace Amanzi

#endif

