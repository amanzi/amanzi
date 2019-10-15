/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon (coonet@ornl.gov)
      Konstantin Lipnikov (lipnikov@lanl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>

template <class Vector, class VectorSpace>
class MatrixBase {
 public:
  // NOTE that a default constructor MUST be allowed, even if it
  // cannot work if actually used.  This is because LinearOperator,
  // when given a Matrix, inherits from that Matrix's class so as to
  // be able to replace that class.  Note that nothing is inherited
  // from the base class, and so no functionality from the
  // default-constructed MatrixBase is ever used.
  MatrixBase();

  // Space for the domain of the operator.
  const VectorSpace& DomainMap() const;

  // Space for the domain of the operator.
  const VectorSpace& RangeMap() const;

  // Apply matrix, b <-- Ax, returns ierr = 0 if success, !0 otherwise
  int Apply(const Vector& x, Vector& b) const;

  // Apply the inverse, x <-- A^-1 b, returns ierr >= 0 if success, <0 otherwise
  int ApplyInverse(const Vector& b, Vector& x) const;
};
