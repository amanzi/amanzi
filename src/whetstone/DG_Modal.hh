/*
  WhetStone, Version 2.2
  Release name: naka-to.

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Discontinuous Galerkin modal method. Efficient implementation
  requires to cache various data for all mesh cells.
*/

#ifndef AMANZI_WHETSTONE_DG_MODAL_HH_
#define AMANZI_WHETSTONE_DG_MODAL_HH_

#include <memory>

#include "Teuchos_RCP.hpp"

// Amanzi
#include "MeshLight.hh"
#include "Point.hh"

// WhetStone
#include "Basis.hh"
#include "BilinearForm.hh"
#include "BilinearFormFactory.hh"
#include "DenseMatrix.hh"
#include "DenseVector.hh"
#include "MatrixObjects.hh"
#include "NumericalIntegration.hh"
#include "Polynomial.hh"
#include "PolynomialOnMesh.hh"
#include "Tensor.hh"
#include "VectorObjects.hh"
#include "WhetStoneDefs.hh"
#include "WhetStoneFunction.hh"

namespace Amanzi {
namespace WhetStone {

class Polynomial;

class DG_Modal : public BilinearForm {
 public:
  DG_Modal(const Teuchos::ParameterList& plist,
           const Teuchos::RCP<const AmanziMesh::MeshLight>& mesh);

  // basic member functions
  // -- schema
  virtual std::vector<SchemaItem> schema() const override;

  // -- mass matrices
  virtual int MassMatrix(int c, const Tensor& K, DenseMatrix& M) override;
  virtual int MassMatrix(int c, const Polynomial& K, DenseMatrix& M) override;
  int MassMatrix(int c, const Tensor& K, PolynomialOnMesh& integrals, DenseMatrix& M);

  // -- stiffness matrices. General coefficient requires to specify the quadrature order
  virtual int StiffnessMatrix(int c, const Tensor& K, DenseMatrix& A) override;
  virtual int StiffnessMatrix(int c, const WhetStoneFunction* K, DenseMatrix& A) override;
  virtual int StiffnessMatrix(int c, const MatrixPolynomial& K, DenseMatrix& A) override;

  // -- advection matrices
  virtual int AdvectionMatrix(int c, const VectorPolynomial& uc,
                              DenseMatrix& A, bool grad_on_test) override;

  // -- flux matrices
  //    returns point flux value (u.n) in the last parameter
  int FluxMatrix(int f, const Polynomial& uf, DenseMatrix& A, bool upwind, bool jump_on_test, double* flux);
  int FluxMatrixRusanov(int f, const VectorPolynomial& uc1, const VectorPolynomial& uc2,
                        const Polynomial& uf, DenseMatrix& A);
  int FluxMatrixGaussPoints(int f, const Polynomial& uf, DenseMatrix& A, bool upwind, bool jump_on_test);

  // -- interface matrices: jumps and penalty
  template<typename Coef, typename std::enable_if<!std::is_pointer<Coef>::value>::type* = nullptr>
  int FaceMatrixJump(int f, const Coef& K1, const Coef& K2, DenseMatrix& A);
  int FaceMatrixJump(int f, const WhetStoneFunction* K1, const WhetStoneFunction* K2, DenseMatrix& A);

  int FaceMatrixPenalty(int f, double Kf, DenseMatrix& A);

  // miscalleneous
  // -- order of polynomials in each cell
  void set_order(int order) { order_ = order; }
  int get_order() const { return order_; }

  // -- access
  const Basis<AmanziMesh::MeshLight>& cell_basis(int c) const { return *basis_[c]; }
  Polynomial& monomial_integrals(int c) { return monomial_integrals_[c]; }

 private:
  int numi_order_;
  NumericalIntegration<AmanziMesh::MeshLight> numi_;

  std::vector<Polynomial> monomial_integrals_;  // integrals of non-normalized monomials
  std::vector<std::shared_ptr<Basis<AmanziMesh::MeshLight> > > basis_;

  static RegisteredFactory<DG_Modal> factory_;
};


/* *****************************************************************
* Jump matrix for Taylor basis using tensors:
*   \Int_f ( {K \grad \rho} [\psi] ) dS
****************************************************************** */
template<typename Coef, typename std::enable_if<!std::is_pointer<Coef>::value>::type*>
int DG_Modal::FaceMatrixJump(int f, const Coef& K1, const Coef& K2, DenseMatrix& A)
{
  AmanziMesh::Entity_ID_List cells;
  mesh_->face_get_cells(f, Parallel_type::ALL, &cells);
  int ncells = cells.size();

  Polynomial poly(d_, order_);
  int size = poly.size();

  int nrows = ncells * size;
  A.Reshape(nrows, nrows);

  // Calculate integrals needed for scaling
  int c1 = cells[0];
  int c2 = (ncells > 1) ? cells[1] : -1;

  // Calculate co-normals
  int dir;
  AmanziGeometry::Point normal = mesh_->face_normal(f, false, c1, &dir);

  normal /= norm(normal);
  auto conormal1 = K1 * normal;
  auto conormal2 = K2 * normal;

  // integrate traces of polynomials on face f
  double coef00, coef01, coef10, coef11;
  VectorPolynomial pgrad;
  std::vector<const PolynomialBase*> polys(2);

  for (auto it = poly.begin(); it < poly.end(); ++it) {
    const int* idx0 = it.multi_index();
    int k = PolynomialPosition(d_, idx0);

    Polynomial p0(d_, idx0, 1.0);
    p0.set_origin(mesh_->cell_centroid(c1));

    pgrad = Gradient(p0);
    p0 = pgrad * conormal1;

    for (auto jt = poly.begin(); jt < poly.end(); ++jt) {
      const int* idx1 = jt.multi_index();
      int l = PolynomialPosition(d_, idx1);

      Polynomial q0(d_, idx1, 1.0);
      q0.set_origin(mesh_->cell_centroid(c1));

      polys[0] = &p0;
      polys[1] = &q0;
      coef00 = numi_.IntegratePolynomialsFace(f, polys);

      A(k, l) = coef00 / ncells;

      if (c2 >= 0) {
        Polynomial p1(d_, idx0, 1.0);
        p1.set_origin(mesh_->cell_centroid(c2));

        pgrad = Gradient(p1);
        p1 = pgrad * conormal2;

        Polynomial q1(d_, idx1, 1.0);
        q1.set_origin(mesh_->cell_centroid(c2));

        polys[1] = &q1;
        coef01 = numi_.IntegratePolynomialsFace(f, polys);

        polys[0] = &p1;
        coef11 = numi_.IntegratePolynomialsFace(f, polys);

        polys[1] = &q0;
        coef10 = numi_.IntegratePolynomialsFace(f, polys);

        A(k, size + l) = -coef01 / ncells;
        A(size + k, size + l) = -coef11 / ncells;
        A(size + k, l) = coef10 / ncells;
      }
    }
  }

  if (ncells == 1) {
    basis_[c1]->BilinearFormNaturalToMy(A);
  } else {
    basis_[c1]->BilinearFormNaturalToMy(basis_[c1], basis_[c2], A);
  }

  return 0;
}

}  // namespace WhetStone
}  // namespace Amanzi

#endif

