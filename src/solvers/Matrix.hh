/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
           Konstantin Lipnikov (lipnikov@lanl.gov)
*/
//! Base class for providing Apply() and applyInverse() methods.

/*!

The following interface is implemented by all Operator-like classes, and
enables both the forward and inverse operators to be applid.

Developer notes:

- Apply(X,Y) requires (but may not check) that X is from DomainMap(), and Y
  from RangeMap().  It returns 0 on success and 1 on failure, and may (or may
  not) throw on error.

- set_inverse_parameters() is frequently just a setter, and likely does not
  instantiate any underlying objects.  It may be called on construction of the
  object, if an object ParameterList is supplied.  It must be called prior to...

- InitializeInverse() instantiates inner objects and may be called once
  symbolic structure is known (e.g. all maps, local ops, block operators, etc
  are set, and SymbolicAssemble() is now valid).  It must be called prior to
  ComputeInverse(), but may be called lazily by ComputeInverse().

- ComputeInverse() does numerical work on the values of the operator
  (e.g. factorizations, etc).  It should be called each time the operator's
  values have changed, and may be called lazily by applyInverse().

- applyInverse(X,Y) requires (but may not check) that X is from RangeMap() and
  Y from DomainMap().  It returns 0 on success and something else on failure,
  for ALL METHODS!

- returned_code() and returned_code_string() may be used to get the actual
  value returned by the inner method implementing applyInverse(), and its
  transliteration.  Frequently, for iterative methods, this is postive and
  equal to the number of iterations if converged, and negative if failed.  For
  direct methods and preconditioners, it may be set by a TPL, or may simply be
  0 on success and 1 on failure, or 0 on success and a random, unparsed,
  unknown value on failure.

 */

#pragma once

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "dbc.hh"
#include "errors.hh"

namespace Amanzi {

template<class Vector,
         class VectorSpace>
class Matrix {
 public:
  using Vector_t = Vector;
  using VectorSpace_t = VectorSpace;

  virtual ~Matrix() = default;

  virtual int apply(const Vector& x, Vector& y) const = 0;

  virtual void set_inverse_parameters(Teuchos::ParameterList& inv_list) = 0;
  virtual void initializeInverse() = 0;
  virtual void computeInverse() = 0;
  virtual int applyInverse(const Vector& y, Vector& x) const = 0;

  virtual const Teuchos::RCP<const VectorSpace> getDomainMap() const = 0;
  virtual const Teuchos::RCP<const VectorSpace> getRangeMap() const = 0;

  //  virtual double TrueResidual(const Vector& x, const Vector& y) const = 0;

  // control and statistics -- must be valid for both iterative and
  // non-iterative methods, approximate and exact methods.
  virtual double residual() const = 0;
  virtual int num_itrs() const = 0;
  //  virtual void add_criteria(int criteria) = 0;

  virtual int returned_code() const = 0;
  virtual std::string returned_code_string() const = 0;

  virtual std::string name() const = 0;
};

}  // namespace Amanzi




