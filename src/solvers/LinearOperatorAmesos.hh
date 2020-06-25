/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Konstantin Lipnikov (amklinv@sandia.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>

#ifndef AMANZI_AMESOS_OPERATOR_HH_
#define AMANZI_AMESOS_OPERATOR_HH_

#include <cmath>

#include "Amesos.h"
#include "Amesos2.hpp"
#include "Amesos_BaseSolver.h"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "CompositeMatrix.hh"
#include "errors.hh"
#include "VerboseObject.hh"

#include "LinearOperator.hh"
#include "LinearOperatorDefs.hh"
#include "MatrixJF.hh"

namespace Amanzi {
namespace AmanziSolvers {

/* ******************************************************************
 * Auxiliary base class.
 ****************************************************************** */
template <class Matrix, class Vector, class VectorSpace>
class LinearOperatorAmesosBase
  : public LinearOperator<Matrix, Vector, VectorSpace> {
 public:
  LinearOperatorAmesosBase(const Teuchos::RCP<const Matrix>& m,
                           const Teuchos::RCP<const Matrix>& h)
    : LinearOperator<Matrix, Vector, VectorSpace>(m, h){};

  void Init(Teuchos::ParameterList& plist);
  void Init() { LinearOperator<Matrix, Vector, VectorSpace>::Init(); }

  // access
  void add_criteria(int criteria){};
  int num_itrs() { return 0; }
  double residual() { return 0.0; }
  int returned_code() { return 0; }

 protected:
  using LinearOperator<Matrix, Vector, VectorSpace>::name_;

  mutable Teuchos::ParameterList plist_;
  Teuchos::RCP<VerboseObject> vo_;
  mutable bool initialized_ = false;
  mutable int version_;
};


/* ******************************************************************
 * Primary interface to direct solvers in Amesos.
 ****************************************************************** */
template <class Matrix, class Vector, class VectorSpace>
class LinearOperatorAmesos
  : public LinearOperatorAmesosBase<Matrix, Vector, VectorSpace> {
 public:
  LinearOperatorAmesos(const Teuchos::RCP<const Matrix>& m,
                       const Teuchos::RCP<const Matrix>& h)
    : LinearOperatorAmesosBase<Matrix, Vector, VectorSpace>(m, h){};

  int ApplyInverse(const Vector& v, Vector& hv) const;

 private:
  int ApplyInverseVer1_(const Vector& v, Vector& hv) const;
  int ApplyInverseVer2_(const Vector& v, Vector& hv) const;

 protected:
  using LinearOperator<Matrix, Vector, VectorSpace>::m_;
  using LinearOperator<Matrix, Vector, VectorSpace>::name_;
  using LinearOperatorAmesosBase<Matrix, Vector, VectorSpace>::plist_;
  using LinearOperatorAmesosBase<Matrix, Vector, VectorSpace>::initialized_;
  using LinearOperatorAmesosBase<Matrix, Vector, VectorSpace>::version_;
};


/* ******************************************************************
 * Partial specification for MatrixJF<Vector, VectorSpace>
 ****************************************************************** */
template <class Vector, class VectorSpace>
class LinearOperatorAmesos<MatrixJF<Vector, VectorSpace>, Vector, VectorSpace>
  : public LinearOperatorAmesosBase<MatrixJF<Vector, VectorSpace>, Vector,
                                    VectorSpace> {
 public:
  LinearOperatorAmesos(
    const Teuchos::RCP<const MatrixJF<Vector, VectorSpace>>& m,
    const Teuchos::RCP<const MatrixJF<Vector, VectorSpace>>& h)
    : LinearOperatorAmesosBase<MatrixJF<Vector, VectorSpace>, Vector,
                               VectorSpace>(m, h){};

  virtual int ApplyInverse(const Vector& v, Vector& hv) const override
  {
    Errors::Message msg("LinearOperatorAmesos: missing partial specification.");
    Exceptions::amanzi_throw(msg);
    return LIN_SOLVER_AMESOS_SAYS_FAIL;
  }
};


/* ******************************************************************
 * Partial specification for CompositeMatrix
 ****************************************************************** */
template <class Vector, class VectorSpace>
class LinearOperatorAmesos<CompositeMatrix, Vector, VectorSpace>
  : public LinearOperatorAmesosBase<CompositeMatrix, Vector, VectorSpace> {
 public:
  LinearOperatorAmesos(const Teuchos::RCP<const CompositeMatrix>& m,
                       const Teuchos::RCP<const CompositeMatrix>& h)
    : LinearOperatorAmesosBase<CompositeMatrix, Vector, VectorSpace>(m, h){};

  virtual int ApplyInverse(const Vector& v, Vector& hv) const override
  {
    Errors::Message msg("LinearOperatorAmesos: missing partial specification.");
    Exceptions::amanzi_throw(msg);
    return LIN_SOLVER_AMESOS_SAYS_FAIL;
  }
};


/* ******************************************************************
 * Initialization from a parameter list.
 ****************************************************************** */
template <class Matrix, class Vector, class VectorSpace>
void
LinearOperatorAmesosBase<Matrix, Vector, VectorSpace>::Init(
  Teuchos::ParameterList& plist)
{
  plist_ = plist;

  vo_ = Teuchos::rcp(new VerboseObject("Solvers::Amesos", plist));
  name_ = plist.get<std::string>("solver name", "Klu");
  version_ = plist.get<int>("amesos version", 1);
  initialized_ = true;

  if (version_ != 1 && version_ != 2) {
    Errors::Message msg("LinearOperatorAmesos: unsupported version of Amesos.");
    Exceptions::amanzi_throw(msg);
  }
}


/* ******************************************************************
 * Generic implementation
 ****************************************************************** */
template <class Matrix, class Vector, class VectorSpace>
int
LinearOperatorAmesos<Matrix, Vector, VectorSpace>::ApplyInverse(
  const Vector& v, Vector& hv) const
{
  if (!initialized_) {
    Errors::Message msg("LinearOperatorAmesos: has not been initialized.");
    Exceptions::amanzi_throw(msg);
  }

  if (version_ == 1) {
    return ApplyInverseVer1_(v, hv);
  } else if (version_ == 2) {
    return ApplyInverseVer2_(v, hv);
  }

  return LIN_SOLVER_AMESOS_SAYS_SUCCESS;
}


/* ******************************************************************
 * Generic implementation for Amesos1
 ****************************************************************** */
template <class Matrix, class Vector, class VectorSpace>
int
LinearOperatorAmesos<Matrix, Vector, VectorSpace>::ApplyInverseVer1_(
  const Vector& v, Vector& hv) const
{
  Epetra_Vector rhs(m_->A()->RowMap());
  Epetra_Vector sol(rhs);

  m_->CopyVectorToSuperVector(v, rhs);

  Epetra_LinearProblem problem(
    const_cast<Epetra_CrsMatrix*>(&*m_->A()), &sol, &rhs);

  Amesos factory;
  auto solver = factory.Create(name_, problem);
  if (solver == NULL) {
    Errors::Message msg;
    msg << "LinearOperatorAmesos: solver \"" << name_ << "\" is not available";
    Exceptions::amanzi_throw(msg);
  }

  solver->SetParameters(plist_);

  int ierr = solver->SymbolicFactorization();
  if (ierr > 0) return LIN_SOLVER_AMESOS_SYMBOLIC_FAIL;

  ierr = solver->NumericFactorization();
  if (ierr > 0) return LIN_SOLVER_AMESOS_FACTORIZATION_FAIL;

  try {
    solver->Solve();
    delete solver;
  } catch (...) {
    Errors::Message msg("LinearOperatorAmesos: solver failed");
    Exceptions::amanzi_throw(msg);
  }

  m_->CopySuperVectorToVector(sol, hv);

  return LIN_SOLVER_AMESOS_SAYS_SUCCESS;
}


/* ******************************************************************
 * Generic implementation for Amesos2
 ****************************************************************** */
template <class Matrix, class Vector, class VectorSpace>
int
LinearOperatorAmesos<Matrix, Vector, VectorSpace>::ApplyInverseVer2_(
  const Vector& v, Vector& hv) const
{
  Epetra_Vector rhs(m_->A()->RowMap());
  Epetra_Vector sol(rhs);

  m_->CopyVectorToSuperVector(v, rhs);

  auto solver = Amesos2::create<Epetra_CrsMatrix, Epetra_MultiVector>(
    name_, m_->A(), Teuchos::rcpFromRef(sol), Teuchos::rcpFromRef(rhs));

  try {
    solver->solve();
  } catch (...) {
    Errors::Message msg("LinearOperatorAmesos: solver failed");
    Exceptions::amanzi_throw(msg);
  }

  m_->CopySuperVectorToVector(sol, hv);

  return LIN_SOLVER_AMESOS_SAYS_SUCCESS;
}

} // namespace AmanziSolvers
} // namespace Amanzi

#endif
