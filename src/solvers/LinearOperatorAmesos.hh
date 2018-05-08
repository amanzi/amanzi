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
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "errors.hh"
#include "VerboseObject.hh"

#include "DenseMatrix.hh"
#include "LinearOperator.hh"
#include "LinearOperatorDefs.hh"

namespace Amanzi {
namespace AmanziSolvers {

template<class Matrix, class Vector, class VectorSpace>
class LinearOperatorAmesos : public LinearOperator<Matrix, Vector, VectorSpace>
{
 public:
  LinearOperatorAmesos(const Teuchos::RCP<const Matrix>& m,
                       const Teuchos::RCP<const Matrix>& h) :
      LinearOperator<Matrix, Vector, VectorSpace>(m, h),
      initialized_(false) {};

  void Init(Teuchos::ParameterList& plist);
  void Init() { LinearOperator<Matrix, Vector, VectorSpace>::Init(); }

  int ApplyInverse(const Vector& v, Vector& hv) const;

  // access
  void add_criteria(int criteria) {};
  int num_itrs() { return 0; }
  double residual() { return 0.0; }
  int returned_code() { return 0; }

 private:
  using LinearOperator<Matrix, Vector, VectorSpace>::m_;
  using LinearOperator<Matrix, Vector, VectorSpace>::name_;

  Teuchos::RCP<VerboseObject> vo_;
  mutable bool initialized_;
};


/* ******************************************************************
 * Initialization from a parameter list.
 ****************************************************************** */
template<class Matrix, class Vector, class VectorSpace>
void LinearOperatorAmesos<Matrix, Vector, VectorSpace>::Init(Teuchos::ParameterList& plist)
{
  vo_ = Teuchos::rcp(new VerboseObject("Solvers::Amesos", plist));
  name_ = plist.get<std::string>("solver name", "Amesos_Klu");
  initialized_ = true;
}


template<class Matrix, class Vector, class VectorSpace>
int LinearOperatorAmesos<Matrix, Vector, VectorSpace>::ApplyInverse(const Vector& v, Vector& hv) const
{
  if (!initialized_) {
    Errors::Message msg("LinearOperatorAmesos: has not been initialized.");
    Exceptions::amanzi_throw(msg);
  }

  Vector rhs(v);
  Epetra_LinearProblem problem(&*m_->A(), &hv, &rhs);

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

}  // namespace AmanziSolvers
}  // namespace Amanzi

#endif
