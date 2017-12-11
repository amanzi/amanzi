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

#include "Mesh.hh"
#include "Point.hh"

#include "DenseMatrix.hh"
#include "InnerProductL2.hh"
#include "InnerProductH1.hh"
#include "Tensor.hh"
#include "VectorPolynomial.hh"

namespace Amanzi {
namespace WhetStone {

class Polynomial;

class BilinearForm : public virtual InnerProductL2,
                     public virtual InnerProductH1 {

 public:
  explicit BilinearForm() {};
  ~BilinearForm() {};

  // additional members
  // -- low-order schemes require constant vector/tensor coefficients
  //    we also specify function to which gradient operator is applied
  virtual int AdvectionMatrix(int c, const AmanziGeometry::Point v, DenseMatrix& A, bool grad_on_test) = 0;
  virtual int DivergenceMatrix(int c, DenseMatrix& A) = 0;

  // -- high-order schemes may require polynomial coefficients
  //    We use ending "Poly" to simplify logic of derived classes (no keyword "using")
  virtual int MassMatrixPoly(int c, const Polynomial& K, DenseMatrix& M) = 0; 
  virtual int StiffnessMatrixPoly(int c, const Polynomial& K, DenseMatrix& A) = 0;
  virtual int AdvectionMatrixPoly(int c, const VectorPolynomial& v, DenseMatrix& A, bool grad_on_test) = 0;
};

}  // namespace WhetStone
}  // namespace Amanzi

#endif

