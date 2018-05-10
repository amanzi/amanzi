/*
  Solvers

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (amklinv@sandia.gov)

  Direct solvers via Trilinos.
*/

#ifndef  AMANZI_AMESOS_OPERATOR_HH_
#define  AMANZI_AMESOS_OPERATOR_HH_

#include <cmath>

#include "Amesos.h"
#include "Amesos_BaseSolver.h"
#include "Epetra_CrsMatrix.h"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "errors.hh"
#include "CompositeVector.hh"
#include "OperatorUtils.hh"
#include "VerboseObject.hh"

#include "DenseMatrix.hh"
#include "LinearOperator.hh"
#include "LinearOperatorDefs.hh"
#include "MatrixJF.hh"

namespace Amanzi {
namespace AmanziSolvers {

/* ******************************************************************
* Auxiliary base class.
****************************************************************** */
template<class Matrix, class Vector, class VectorSpace>
class LinearOperatorAmesosBase : public LinearOperator<Matrix, Vector, VectorSpace> {
 public:
  LinearOperatorAmesosBase(const Teuchos::RCP<const Matrix>& m,
                           const Teuchos::RCP<const Matrix>& h) :
      LinearOperator<Matrix, Vector, VectorSpace>(m, h) {};

  void Init(Teuchos::ParameterList& plist);
  void Init() { LinearOperator<Matrix, Vector, VectorSpace>::Init(); }

  // virtual int ApplyInverse(const Vector& v, Vector& hv) const;

  // access
  void add_criteria(int criteria) {};
  int num_itrs() { return 0; }
  double residual() { return 0.0; }
  int returned_code() { return 0; }

 protected:
  using LinearOperator<Matrix, Vector, VectorSpace>::m_;
  using LinearOperator<Matrix, Vector, VectorSpace>::name_;

  Teuchos::RCP<VerboseObject> vo_;
  mutable bool initialized_ = false;
};


/* ******************************************************************
* Initialization from a parameter list.
****************************************************************** */
template<class Matrix, class Vector, class VectorSpace>
void LinearOperatorAmesosBase<Matrix, Vector, VectorSpace>::Init(Teuchos::ParameterList& plist)
{
  vo_ = Teuchos::rcp(new VerboseObject("Solvers::Amesos", plist));
  name_ = plist.get<std::string>("solver name", "Amesos_Klu");
  initialized_ = true;
}


/* ******************************************************************
* Primary interface to dierct solvers in Amesos.
****************************************************************** */
template<class Matrix, class Vector, class VectorSpace>
class LinearOperatorAmesos : public LinearOperatorAmesosBase<Matrix, Vector, VectorSpace>
{
 public:
  LinearOperatorAmesos(const Teuchos::RCP<const Matrix>& m,
                       const Teuchos::RCP<const Matrix>& h) :
      LinearOperatorAmesosBase<Matrix, Vector, VectorSpace>(m, h) {};

  // default implementation is empty for the moment
  virtual int ApplyInverse(const Vector& v, Vector& hv) const { 
    Errors::Message msg("LinearOperatorAmesos: missing partial specification.");
    Exceptions::amanzi_throw(msg);
    return LIN_SOLVER_AMESOS_SAYS_FAIL; 
  }
};


/* ******************************************************************
* Partial specification for CompositeVector and CompositeVectorSpace
****************************************************************** */
template<class Matrix> 
class LinearOperatorAmesos<Matrix, CompositeVector, CompositeVectorSpace> :
    public LinearOperatorAmesosBase<Matrix, CompositeVector, CompositeVectorSpace> {
 public:
  LinearOperatorAmesos(const Teuchos::RCP<const Matrix>& m,
                       const Teuchos::RCP<const Matrix>& h) :
      LinearOperatorAmesosBase<Matrix, CompositeVector, CompositeVectorSpace>(m, h) {};

  int ApplyInverse(const CompositeVector& v, CompositeVector& hv) const {
    if (!initialized_) {
      Errors::Message msg("LinearOperatorAmesos: has not been initialized.");
      Exceptions::amanzi_throw(msg);
    }

    Epetra_Vector rhs(m_->A()->RowMap());
    Epetra_Vector sol(rhs);

    Operators::CopyCompositeVectorToSuperVector(*m_->smap(), v, rhs, m_->schema_col());
    Epetra_LinearProblem problem(const_cast<Epetra_CrsMatrix*>(&*m_->A()), &sol, &rhs);

    Amesos factory;
    auto solver = factory.Create(name_, problem);
    if (solver == NULL) {
      Errors::Message msg("LinearOperatorAmesos: specified solver is not available");
      Exceptions::amanzi_throw(msg);
    }

    int ierr = solver->SymbolicFactorization();
    if (ierr > 0) return LIN_SOLVER_AMESOS_SYMBOLIC_FAIL;

    ierr = solver->NumericFactorization();
    if (ierr > 0) return LIN_SOLVER_AMESOS_FACTORIZATION_FAIL;

    try {
      solver->Solve();
    } catch(...) {
      Errors::Message msg("LinearOperatorAmesos: solver failed");
      Exceptions::amanzi_throw(msg);
    }

    Operators::CopySuperVectorToCompositeVector(*m_->smap(), sol, hv, m_->schema_col());

    delete solver;
    return LIN_SOLVER_AMESOS_SAYS_SUCCESS;
  }

 protected:
  using LinearOperatorAmesosBase<Matrix, CompositeVector, CompositeVectorSpace>::m_;
  using LinearOperatorAmesosBase<Matrix, CompositeVector, CompositeVectorSpace>::name_;
  using LinearOperatorAmesosBase<Matrix, CompositeVector, CompositeVectorSpace>::initialized_;
};


/* ******************************************************************
* Partial specification for Epetra_MultiVector and Epetra_Map
****************************************************************** */
template<class Matrix>
class LinearOperatorAmesos<Matrix, Epetra_MultiVector, Epetra_Map> :
    public LinearOperatorAmesosBase<Matrix, Epetra_MultiVector, Epetra_Map> {
 public:
  LinearOperatorAmesos(const Teuchos::RCP<const Matrix>& m,
                       const Teuchos::RCP<const Matrix>& h) :
      LinearOperatorAmesosBase<Matrix, Epetra_MultiVector, Epetra_Map>(m, h) {};

  int ApplyInverse(const Epetra_MultiVector& v, Epetra_MultiVector& hv) const {
    if (!initialized_) {
      Errors::Message msg("LinearOperatorAmesos: has not been initialized.");
      Exceptions::amanzi_throw(msg);
    }

    Epetra_MultiVector rhs(v);
    Epetra_LinearProblem problem(const_cast<Epetra_CrsMatrix*>(&*m_->A()), &hv, &rhs);

    Amesos factory;
    auto solver = factory.Create(name_, problem);
    if (solver == NULL) {
      Errors::Message msg("LinearOperatorAmesos: specified solver is not available");
      Exceptions::amanzi_throw(msg);
    }

    int ierr = solver->SymbolicFactorization();
    if (ierr > 0) return LIN_SOLVER_AMESOS_SYMBOLIC_FAIL;

    ierr = solver->NumericFactorization();
    if (ierr > 0) return LIN_SOLVER_AMESOS_FACTORIZATION_FAIL;

    solver->Solve();

    delete solver;
    return LIN_SOLVER_AMESOS_SAYS_SUCCESS;
  }

 protected:
  using LinearOperatorAmesosBase<Matrix, Epetra_MultiVector, Epetra_Map>::m_;
  using LinearOperatorAmesosBase<Matrix, Epetra_MultiVector, Epetra_Map>::name_;
  using LinearOperatorAmesosBase<Matrix, Epetra_MultiVector, Epetra_Map>::initialized_;
};


/* ******************************************************************
* Full specification for MatrixJF, CompositeVector, CompositeVectorSpace
****************************************************************** */
typedef MatrixJF<CompositeVector, CompositeVectorSpace> MatrixJFcv;
template<>
class LinearOperatorAmesos<MatrixJFcv, CompositeVector, CompositeVectorSpace> :
    public LinearOperatorAmesosBase<MatrixJFcv, CompositeVector, CompositeVectorSpace> {
 public:
  LinearOperatorAmesos(const Teuchos::RCP<const MatrixJFcv>& m,
                       const Teuchos::RCP<const MatrixJFcv>& h) :
      LinearOperatorAmesosBase<MatrixJFcv, CompositeVector, CompositeVectorSpace>(m, h) {};

  int ApplyInverse(const CompositeVector& v, CompositeVector& hv) const {
    Errors::Message msg("LinearOperatorAmesos: missing partial specification.");
    Exceptions::amanzi_throw(msg);
    return LIN_SOLVER_AMESOS_SAYS_FAIL; 
  }

 protected:
  using LinearOperatorAmesosBase<MatrixJFcv, CompositeVector, CompositeVectorSpace>::m_;
  using LinearOperatorAmesosBase<MatrixJFcv, CompositeVector, CompositeVectorSpace>::name_;
  using LinearOperatorAmesosBase<MatrixJFcv, CompositeVector, CompositeVectorSpace>::initialized_;
};
}  // namespace AmanziSolvers
}  // namespace Amanzi

#endif
