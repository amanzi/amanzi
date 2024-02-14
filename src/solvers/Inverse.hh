/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
           Konstantin Lipnikov (lipnikov@lanl.gov)
*/

//! Base class for providing applyInverse() methods on operators.
/*!

Matrix provides existing Operators with an inverse.  Note this may be
iterative or non-iterative, assembled or non-assembled, approximate or exact to
machine precision.

.. _inverse-typed-spec:
.. admonition:: inverse-typed-spec

   * `"iterative method`" ``[string]`` **optional**
   * `"direct method`" ``[string]`` **optional**
   * `"preconditioning method`" ``[string]`` **optional**

DOCUMENT ME!

*/

/*
Developer notes:
----------------

Matrix takes a four-stage approach, following the more modern Trilinos
packages (init(), update(), compute(), apply()):

- set_inverse_parameters() processes the ParameterList, parsing options.

- InitializeInverse() implies that the symbolic structure is now known.  Changes
to symbolic structure require calling InitializeInverse() again.  All work that
can leverage this, e.g. allocation of work space, etc, can now be done.

- ComputeInverse() requires that values in the operator have now been set.
  Work such as calculating L and U, etc, can now be done.

- applyInverse() accepts vectors and applies the inverse.  It returns 0 on
  success and 1 on failure.

Note that any stage may be called without invalidating any stage before it, but
necessitates calling all stages after it.

*/
#pragma once

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "dbc.hh"
#include "errors.hh"

#include "Matrix.hh"
#include "InverseDefs.hh"

namespace Amanzi {
namespace AmanziSolvers {

template <class Operator,
          class Preconditioner = Operator,
          class Vector = typename Operator::Vector_t,
          class VectorSpace = typename Vector::VectorSpace_t>
class Inverse : public Matrix<Vector, VectorSpace> {
 public:
  Inverse() = default;
  Inverse(const Inverse& other) = delete;
  virtual ~Inverse() = default;

  virtual void set_matrices(const Teuchos::RCP<Operator>& m, const Teuchos::RCP<Preconditioner>& h)
  {
    m_ = m;
    h_ = h;
  }

  template <typename T = Preconditioner>
  typename std::enable_if<std::is_same<Operator, T>::value>::type
  set_matrix(const Teuchos::RCP<Operator>& m)
  {
    set_matrices(m, m);
  }

  virtual const Teuchos::RCP<const VectorSpace> getDomainMap() const { return m_->getDomainMap(); }
  virtual const Teuchos::RCP<const VectorSpace> getRangeMap() const { return m_->getRangeMap(); }

  int apply(const Vector& x, Vector& y) const
  {
    m_->apply(x, y);
    return 0;
  }
  virtual void set_inverse_parameters(Teuchos::ParameterList& plist) = 0;
  virtual void initializeInverse() = 0;
  virtual void computeInverse() = 0;

  virtual int applyInverse(const Vector& X, Vector& Y) const = 0;

  // This ensures that anything that is an Inverse can be used as a
  // preconditioner.
  int ApplyInverseUserSupplied(const Vector& X, Vector& Y) const { return applyInverse(X, Y); }

  double TrueResidual(const Vector& x, const Vector& y) const
  {
    Vector r(y);
    m_->apply(x, r); // r = y - M * x
    r.update(1.0, y, -1.0);

    return r.norm2();
  }

  // control and statistics -- must be valid for both iterative and
  // non-iterative methods, approximate and exact methods.
  virtual double residual() const { return 0.; }
  virtual int num_itrs() const { return 0; }
  virtual void add_criteria(int criteria) {}

  virtual int returned_code() const = 0;
  virtual std::string returned_code_string() const = 0;

  std::string name() const { return name_; }
  void set_name(const std::string& name) { name_ = name; }

 protected:
  Teuchos::RCP<Operator> m_;
  Teuchos::RCP<Preconditioner> h_;
  std::string name_;
};

} // namespace AmanziSolvers
} // namespace Amanzi
