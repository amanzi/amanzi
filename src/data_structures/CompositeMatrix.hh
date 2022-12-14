/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/* -------------------------------------------------------------------------

ATS

Interface for a Matrix that acts on CompositeVector.
------------------------------------------------------------------------- */

#ifndef AMANZI_COMPOSITEMATRIX_HH_
#define AMANZI_COMPOSITEMATRIX_HH_

#include "Teuchos_RCP.hpp"

namespace Amanzi {

class CompositeVector;
class CompositeVectorSpace;

class CompositeMatrix {
 public:
  virtual ~CompositeMatrix() = default;

  // Vector space of the Matrix's domain.
  virtual const CompositeVectorSpace& DomainMap() const = 0;

  // Vector space of the Matrix's range.
  virtual const CompositeVectorSpace& RangeMap() const = 0;

  // Apply matrix, b <-- Ax, returns ierr
  virtual int Apply(const CompositeVector& x, CompositeVector& b) const = 0;

  // Apply the inverse, x <-- A^-1 b, returns ierr
  virtual int ApplyInverse(const CompositeVector& b, CompositeVector& x) const = 0;
};

} // namespace Amanzi

#endif
