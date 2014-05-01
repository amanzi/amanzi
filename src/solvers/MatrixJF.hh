/*
  Authors: Ethan Coon (ecoon@lanl.gov)

  This decorator class wraps a nonlinear SolverFnBase as a Matrix
  class that approximates the linear problem at a point by doing a
  standard finite difference (i.e. Knoll and Keyes '04) to approximate
  the action of the Jacobian.  When used with, for instance, NKA, this
  residual evaluation reduces to GMRES on the linear system, giving as
  a solution the standard JFNK correction to the nonlinear problem.

*/


#ifndef AMANZI_JF_MATRIX_HH_
#define AMANZI_JF_MATRIX_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "SolverFnBase.hh"

namespace Amanzi {
namespace AmanziSolvers {

template<class Vector, class VectorSpace>
class MatrixJF {

 public:

  MatrixJF(Teuchos::ParameterList& plist,
	   const Teuchos::RCP<SolverFnBase<Vector> > fn,
	   const VectorSpace& map);

  // Space for the domain of the operator.
  const VectorSpace& DomainMap() const { return *map_; }

  // Space for the domain of the operator.
  const VectorSpace& RangeMap() const { return *map_; }

  // Apply matrix, b <-- Ax, returns ierr = 0 if success, !0 otherwise
  int Apply(const Vector& x, Vector& b) const;

  // Apply the inverse, x <-- A^-1 b, returns ierr = 0 if success, !0 otherwise
  int ApplyInverse(const Vector& b, Vector& x) const;

protected:
  Teuchos::RCP<const VectorSpace> map_;
  Teuchos::RCP<SolverFnBase<Vector> > fn_;
  Teuchos::ParameterList plilst_;
};



} // namespace
} // namespace

#endif

