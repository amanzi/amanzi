/*
  Operators

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  The helper advection-based base class for various remap methods.
*/

// TPLs
#include "Epetra_Vector.h"

// Amanzi
#include "CompositeVector.hh"
#include "TreeVector.hh"

// Amanzi::Operators
#include "RemapDG.hh"

namespace Amanzi {
namespace Operators {

/* *****************************************************************
* Specialization for CompositeVector: Functional evaluation
***************************************************************** */
template<>
void RemapDG<CompositeVector>::FunctionalTimeDerivative(
    double t, const CompositeVector& u, CompositeVector& f)
{
  // -- populate operators
  op_adv_->Setup(velc_, false);
  op_adv_->UpdateMatrices(t);

  op_flux_->Setup(velf_.ptr(), false);
  op_flux_->UpdateMatrices(t);
  op_flux_->ApplyBCs(true, true, true);

  // -- calculate right-hand_side
  op_flux_->global_operator()->Apply(*field_, f);

  nfun_++;
}


/* *****************************************************************
* Limiting the non-conservative field at time t
***************************************************************** */
template<>
void RemapDG<CompositeVector>::ModifySolution(double t, CompositeVector& u)
{
  // populate operators
  op_reac_->Setup(det_, false);
  op_reac_->UpdateMatrices(t);

  // solve the problem with the mass matrix
  auto& matrices = op_reac_->local_op()->matrices;
  for (int n = 0; n < matrices.size(); ++n) matrices[n].InverseSPD();
  op_reac_->global_operator()->Apply(u, *field_);

  // limit non-conservative field and update the conservative field
  if (is_limiter_) {
    // -- save original field and limit it
    auto& field_c = *field_->ViewComponent("cell");
    auto orig_c = field_c;

    ApplyLimiter(t, *field_);

    // -- recover original mass matrices FIXME (lipnikov@lanl.gov)
    for (int n = 0; n < matrices.size(); ++n) matrices[n].InverseSPD();

    // -- shift mean values
    auto& climiter = *limiter_->limiter();
    auto& u_c = *u.ViewComponent("cell");
    int nk = u_c.NumVectors();

    for (int c = 0; c < ncells_owned_; ++c) {
      double a = climiter[c];
      if (a < 1.0) {
        double mass(0.0);
        for (int i = 0; i < nk; ++i) {
          mass += matrices[c](i, 0) * orig_c[i][c];
        }

        field_c[0][c] = a * orig_c[0][c] + (1.0 - a) * mass / matrices[c](0, 0);
      }
    }

    // -- update conservative field
    op_reac_->global_operator()->Apply(*field_, u);
  }
}


/* *****************************************************************
* Change between conservative and non-conservative variables.
***************************************************************** */
template<>
void RemapDG<CompositeVector>::NonConservativeToConservative(
    double t, const CompositeVector& u, CompositeVector& v)
{
  op_reac_->Setup(det_, false);
  op_reac_->UpdateMatrices(t);
  op_reac_->global_operator()->Apply(u, v);
}


template<>
void RemapDG<CompositeVector>::ConservativeToNonConservative(
    double t, const CompositeVector& u, CompositeVector& v)
{
  op_reac_->Setup(det_, false);
  op_reac_->UpdateMatrices(t);

  auto& matrices = op_reac_->local_op()->matrices;
  for (int n = 0; n < matrices.size(); ++n) matrices[n].InverseSPD();

  op_reac_->global_operator()->Apply(u, v);
}


/* *****************************************************************
* Specialization for TreeVector: functional evaluation
***************************************************************** */
template<>
void RemapDG<TreeVector>::FunctionalTimeDerivative(
    double t, const TreeVector& u, TreeVector& f)
{
  // mass conservation equation
  // -- populate operators
  op_adv_->Setup(velc_, false);
  op_adv_->UpdateMatrices(t);

  op_flux_->Setup(velf_.ptr(), false);
  op_flux_->UpdateMatrices(t);
  op_flux_->ApplyBCs(true, true, true);

  // -- calculate right-hand_side
  op_flux_->global_operator()->Apply(*field_, *f.SubVector(0)->Data());

  // volume conservation equation
  auto ones(*field_);
  ones.PutScalar(0.0);
  (*ones.ViewComponent("cell"))(0)->PutScalar(1.0);

  auto tmp = *f.SubVector(1)->Data();
  op_flux_->global_operator()->Apply(ones, tmp);

  // op_reac_->Setup(det_, false);
  // op_reac_->UpdateMatrices(0.0);
  op_reac_->Setup(Teuchos::null);
  op_reac_->UpdateMatrices(Teuchos::null, Teuchos::null);

  auto& matrices = op_reac_->local_op()->matrices;
  for (int n = 0; n < matrices.size(); ++n) matrices[n].InverseSPD();
  op_reac_->global_operator()->Apply(tmp, *f.SubVector(1)->Data());

  nfun_++;
}


/* *****************************************************************
* Limiting the non-conservative field at time t
***************************************************************** */
template<>
void RemapDG<TreeVector>::ModifySolution(double t, TreeVector& u)
{
  // populate operators
  auto detc = *u.SubVector(1)->Data()->ViewComponent("cell");
  int nk = detc.NumVectors();
  WhetStone::DenseVector data(nk);

  for (int c = 0; c < ncells_owned_; ++c) {
    for (int i = 0; i < nk; ++i) data(i) = detc[i][c];
    (*jac_)[c] = dg_->cell_basis(c).CalculatePolynomial(mesh0_, c, order_, data);
  }

  // discrete volume conservation law: new approach
  op_reac_->Setup(jac_, false);
  op_reac_->UpdateMatrices(Teuchos::null, Teuchos::null);

  // discrete volume conservation law: old approach
  // op_reac_->Setup(det_, false);
  // op_reac_->UpdateMatrices(t);

  // solve the problem with the mass matrix
  auto& matrices = op_reac_->local_op()->matrices;
  for (int n = 0; n < matrices.size(); ++n) matrices[n].InverseSPD();
  op_reac_->global_operator()->Apply(*u.SubVector(0)->Data(), *field_);

  // limit non-conservative field and update the conservative field
  if (is_limiter_) {
    // -- save original field and limit it
    auto& field_c = *field_->ViewComponent("cell");
    auto orig_c = field_c;

    ApplyLimiter(t, *field_);

    // -- recover original mass matrices FIXME (lipnikov@lanl.gov)
    for (int n = 0; n < matrices.size(); ++n) matrices[n].InverseSPD();

    // -- shift mean values
    auto& climiter = *limiter_->limiter();
    auto& u_c = *u.SubVector(0)->Data()->ViewComponent("cell");
    int mk = u_c.NumVectors();

    for (int c = 0; c < ncells_owned_; ++c) {
      double a = climiter[c];
      if (a < 1.0) {
        double mass(0.0);
        for (int i = 0; i < mk; ++i) {
          mass += matrices[c](i, 0) * orig_c[i][c];
        }

        field_c[0][c] = a * orig_c[0][c] + (1.0 - a) * mass / matrices[c](0, 0);
      }
    }

    // -- update conservative field
    op_reac_->global_operator()->Apply(*field_, *u.SubVector(0)->Data());
  }
}


/* *****************************************************************
* Change between conservative and non-conservative variables.
***************************************************************** */
template<>
void RemapDG<TreeVector>::NonConservativeToConservative(
    double t, const TreeVector& u, TreeVector& v)
{
  // create a polynomial for determinant of Jacobian
  auto detc = *u.SubVector(1)->Data()->ViewComponent("cell");
  int nk = detc.NumVectors();
  WhetStone::DenseVector data(nk);

  for (int c = 0; c < ncells_owned_; ++c) {
    for (int i = 0; i < nk; ++i) data(i) = detc[i][c];
    (*jac_)[c] = dg_->cell_basis(c).CalculatePolynomial(mesh0_, c, order_, data);
  }

  op_reac_->Setup(jac_, false);
  op_reac_->UpdateMatrices(Teuchos::null, Teuchos::null);

  // conversion is matrix-vector product
  op_reac_->global_operator()->Apply(*u.SubVector(0)->Data(), *v.SubVector(0)->Data());
}


template<>
void RemapDG<TreeVector>::ConservativeToNonConservative(
    double t, const TreeVector& u, TreeVector& v)
{
  // create a polynomial for determinant of Jacobian
  auto detc = *u.SubVector(1)->Data()->ViewComponent("cell");
  int nk = detc.NumVectors();
  WhetStone::DenseVector data(nk);

  for (int c = 0; c < ncells_owned_; ++c) {
    for (int i = 0; i < nk; ++i) data(i) = detc[i][c];
    (*jac_)[c] = dg_->cell_basis(c).CalculatePolynomial(mesh0_, c, order_, data);
  }

  op_reac_->Setup(jac_, false);
  op_reac_->UpdateMatrices(Teuchos::null, Teuchos::null);

  // conversion is inverse matrix-vector product
  auto& matrices = op_reac_->local_op()->matrices;
  for (int n = 0; n < matrices.size(); ++n) matrices[n].InverseSPD();

  op_reac_->global_operator()->Apply(*u.SubVector(0)->Data(), *v.SubVector(0)->Data());
}

}  // namespace Operators
}  // namespace Amanzi

