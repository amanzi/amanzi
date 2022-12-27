/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  WhetStone, Version 2.2
  Release name: naka-to.

  The base virtual class for of bilinear forms of type
     <coef Op1(u), Op2(v)>
  where Op1 and Op2 are differential operators. The tuple
  (u, v, coef) may support a few bilinear forms appearing in
  applications. Most times u = v.
*/

#ifndef AMANZI_BILINEAR_FORM_HH_
#define AMANZI_BILINEAR_FORM_HH_

#include "Teuchos_RCP.hpp"

#include "errors.hh"
#include "MeshLight.hh"

#include "DenseMatrix.hh"
#include "MatrixObjects.hh"
#include "VectorObjects.hh"
#include "WhetStoneDefs.hh"

namespace Amanzi {
namespace WhetStone {

class Polynomial;

class BilinearForm {
 public:
  BilinearForm(const Teuchos::RCP<const AmanziMesh::MeshLight>& mesh)
    : mesh_(mesh), d_(mesh->space_dimension()), order_(1){};
  virtual ~BilinearForm(){};

  // schema
  virtual std::vector<SchemaItem> schema() const = 0;

  // types of bilinear forms differ by the operators applied to arguments
  // -- mass
  virtual int MassMatrix(int c, const Tensor& T, DenseMatrix& M)
  {
    Errors::Message msg("Mass operator with tensor coefficient is not supported.");
    Exceptions::amanzi_throw(msg);
    return 1;
  }

  virtual int MassMatrix(int c, const Polynomial& K, DenseMatrix& M)
  {
    Errors::Message msg("Mass operator with polynomial coefficient is not supported.");
    Exceptions::amanzi_throw(msg);
    return 1;
  }

  virtual int MassMatrixInverse(int c, const Tensor& T, DenseMatrix& W)
  {
    MassMatrix(c, T, W);
    W.Inverse();
    return 0;
  }

  virtual int MassMatrixInverse(int c, const Polynomial& K, DenseMatrix& W)
  {
    MassMatrix(c, K, W);
    W.Inverse();
    return 0;
  }

  // -- stiffness
  virtual int StiffnessMatrix(int c, const Tensor& T, DenseMatrix& A)
  {
    Errors::Message msg("Stiffness operator with tensor is not supported.");
    Exceptions::amanzi_throw(msg);
    return 1;
  }

  virtual int StiffnessMatrix(int c, const Polynomial& K, DenseMatrix& A)
  {
    Errors::Message msg("Stiffness operator with polynomial coefficient is not supported.");
    Exceptions::amanzi_throw(msg);
    return 1;
  }

  virtual int StiffnessMatrix(int c, const MatrixPolynomial& K, DenseMatrix& A)
  {
    Errors::Message msg("Stiffness operator with matrix polynomial coefficient is not supported.");
    Exceptions::amanzi_throw(msg);
    return 1;
  }

  virtual int StiffnessMatrix(int c, const WhetStoneFunction* K, DenseMatrix& A)
  {
    Errors::Message msg("Stiffness operator with function coefficient is not supported.");
    Exceptions::amanzi_throw(msg);
    return 1;
  }

  // -- advection
  virtual int
  AdvectionMatrix(int c, const AmanziGeometry::Point& v, DenseMatrix& A, bool grad_on_test)
  {
    Errors::Message msg("Advection operator is not supported.");
    Exceptions::amanzi_throw(msg);
    return 1;
  }

  virtual int AdvectionMatrix(int c, const VectorPolynomial& v, DenseMatrix& A, bool grad_on_test)
  {
    Errors::Message msg("Advection operator with vector polynomial coefficient is not supported.");
    Exceptions::amanzi_throw(msg);
    return 1;
  }

  virtual int AdvectionMatrix(int c, const std::vector<AmanziGeometry::Point>& u, DenseMatrix& A)
  {
    Errors::Message msg("Advection operator with virtual nodal function is not supported.");
    Exceptions::amanzi_throw(msg);
    return 1;
  }

  // -- divergence
  virtual int DivergenceMatrix(int c, DenseMatrix& A)
  {
    Errors::Message msg("Divergence operator is not supported.");
    Exceptions::amanzi_throw(msg);
    return 1;
  }

  // Projectors is a special type of bilinear form (P(u), v) = (u, v) for any v
  // where coef = 1. If different operators are supported by the derived class,
  // different projectors can be computed.
  // NOTE: we need both (1) action on a function u and (2) matrix form of the
  //       project; hence, the interface will be extended.
  // -- L2 projectors
  virtual void L2Cell(int c,
                      const std::vector<Polynomial>& ve,
                      const std::vector<Polynomial>& vf,
                      const Polynomial* moments,
                      Polynomial& vc)
  {
    Errors::Message msg("L2 projector is not supported.");
    Exceptions::amanzi_throw(msg);
  }
  virtual void
  L2Face(int f, const std::vector<Polynomial>& ve, const Polynomial* moments, Polynomial& vf)
  {
    Errors::Message msg("L2 face projector is not supported.");
    Exceptions::amanzi_throw(msg);
  }

  virtual void L2Cell(int c, const DenseVector& dofs, Polynomial& vc)
  {
    Errors::Message msg("L2 projector (from DOFs) is not supported.");
    Exceptions::amanzi_throw(msg);
  }

  // -- H1 projectors
  virtual void H1Cell(int c,
                      const std::vector<Polynomial>& ve,
                      const std::vector<Polynomial>& vf,
                      const Polynomial* moments,
                      Polynomial& vc)
  {
    Errors::Message msg("H1 cell projector is not supported.");
    Exceptions::amanzi_throw(msg);
  }
  virtual void
  H1Face(int f, const std::vector<Polynomial>& ve, const Polynomial* moments, Polynomial& vf)
  {
    Errors::Message msg("H1 face projector is not supported.");
    Exceptions::amanzi_throw(msg);
  }

  virtual void H1Cell(int c, const DenseVector& dofs, Polynomial& vc)
  {
    Errors::Message msg("H1 cell projector (from DOFs) is not supported.");
    Exceptions::amanzi_throw(msg);
  }

  // miscalleneous
  int get_order() const { return order_; }
  void set_order(int order) { order_ = order; }

 protected:
  Teuchos::RCP<const AmanziMesh::MeshLight> mesh_;
  int d_, order_;
};

} // namespace WhetStone
} // namespace Amanzi

#endif
