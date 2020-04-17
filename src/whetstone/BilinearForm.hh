/*
  WhetStone, Version 2.2
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
#include "MatrixObjects.hh"
#include "Polynomial.hh"
#include "Tensor.hh"
#include "VectorObjects.hh"
#include "WhetStoneDefs.hh"

namespace Amanzi {
namespace WhetStone {

class Polynomial;

class BilinearForm : public virtual InnerProductL2,
                     public virtual InnerProductH1 {
 public:
  explicit BilinearForm(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh)
    : mesh_(mesh),
      d_(mesh->space_dimension()),
      order_(1) {};
  virtual ~BilinearForm() {};

  // schema
  virtual std::vector<SchemaItem> schema() const = 0;

  // additional members
  // -- low-order schemes require typically constant vector/tensor coefficients
  //    also specify function to which gradient operator is applied
  virtual int AdvectionMatrix(int c, const AmanziGeometry::Point v, DenseMatrix& A, bool grad_on_test) {
    Errors::Message msg("AdvectionMatrix: scalar velocity is not supported.");
    Exceptions::amanzi_throw(msg);
    return 1;
  }
  virtual int AdvectionMatrix(int c, const VectorPolynomial& v, DenseMatrix& A, bool grad_on_test) {
    Errors::Message msg("AdvectionMatrix: polynomial coefficient is not supported.");
    Exceptions::amanzi_throw(msg);
    return 1;
  }

  virtual int DivergenceMatrix(int c, DenseMatrix& A) {
    Errors::Message msg("Function DivergenceMatrix is not supported.");
    Exceptions::amanzi_throw(msg);
    return 1;
  }

  // extend interface for the existing members
  // -- high-order schemes may require polynomial coefficients
  using InnerProductL2::MassMatrix;
  virtual int MassMatrix(int c, const Polynomial& K, DenseMatrix& M) {
    Errors::Message msg("MassMatrix: polynomial coefficient is not supported.");
    Exceptions::amanzi_throw(msg);
    return 1;
  }

  using InnerProductL2::MassMatrixInverse;
  virtual int MassMatrixInverse(int c, const Polynomial& K, DenseMatrix& M) {
    Errors::Message msg("MassMatrixInverse: polynomial coefficient is not supported.");
    Exceptions::amanzi_throw(msg);
    return 1;
  }

  using InnerProductH1::StiffnessMatrix;
  virtual int StiffnessMatrix(int c, const Polynomial& K, DenseMatrix& A) {
    Errors::Message msg("StiffnessMatrix: polynomial coefficient is not supported.");
    Exceptions::amanzi_throw(msg);
    return 1;
  }
  // -- general coefficient. Qudrature rule is provided via the input parameter list
  virtual int StiffnessMatrix(int c, const MatrixPolynomial& K, DenseMatrix& A) {
    Errors::Message msg("StiffnessMatrix: matrix polynomial coefficient is not supported.");
    Exceptions::amanzi_throw(msg);
    return 1;
  }
  // -- general coefficient. Qudrature rule is provided via the input parameter list
  virtual int StiffnessMatrix(int c, const WhetStoneFunction* K, DenseMatrix& A) {
    Errors::Message msg("StiffnessMatrix: general coefficient is not supported.");
    Exceptions::amanzi_throw(msg);
    return 1;
  }

  // Projectors facilitate construction of bilinear forms
  // -- L2 projectors
  virtual void L2Cell(int c, const std::vector<Polynomial>& ve,
                      const std::vector<Polynomial>& vf,
                      const Polynomial* moments, Polynomial& vc) {
    Errors::Message msg("L2 projector is not implemented for this scheme.");
    Exceptions::amanzi_throw(msg);
  }
  virtual void L2Face(int f, const std::vector<Polynomial>& ve,
                      const Polynomial* moments, Polynomial& vf) {
    Errors::Message msg("L2 face projector is not implemented for this scheme.");
    Exceptions::amanzi_throw(msg);
  }

  virtual void L2Cell(int c, const DenseVector<>& dofs, Polynomial& vc) {
    Errors::Message msg("L2 projector (from DOFs) is not implemented for this scheme.");
    Exceptions::amanzi_throw(msg);
  }

  // -- H1 projectors
  virtual void H1Cell(int c, const std::vector<Polynomial>& ve,
                      const std::vector<Polynomial>& vf,
                      const Polynomial* moments, Polynomial& vc) {
    Errors::Message msg("H1 cell projector is not implemented for this scheme.");
    Exceptions::amanzi_throw(msg);
  }
  virtual void H1Face(int f, const std::vector<Polynomial>& ve,
                      const Polynomial* moments, Polynomial& vf) {
    Errors::Message msg("H1 face projector is not implemented for this scheme.");
    Exceptions::amanzi_throw(msg);
  }

  virtual void H1Cell(int c, const DenseVector<>& dofs, Polynomial& vc) {
    Errors::Message msg("H1 cell projector (from DOFs) is not implemented for this scheme.");
    Exceptions::amanzi_throw(msg);
  }

  // miscalleneous
  int order() const { return order_; }
  void set_order(int order) { order_ = order; }

 protected:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  int d_, order_;
};

}  // namespace WhetStone
}  // namespace Amanzi

#endif
