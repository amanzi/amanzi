/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------

ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon (ecoon@lanl.gov)

Interface for a Matrix that acts on TreeVector.
------------------------------------------------------------------------- */

#ifndef AMANZI_TREEMATRIX_HH_
#define AMANZI_TREEMATRIX_HH_

#include "Teuchos_RCP.hpp"

namespace Amanzi {

class TreeVector;
class TreeVectorSpace;

class TreeMatrix {
 public:
  // Vector space of the Matrix's domain.
  virtual const TreeVectorSpace& DomainMap() const = 0;

  // Vector space of the Matrix's range.
  virtual const TreeVectorSpace& RangeMap() const = 0;

  // Apply matrix, b <-- Ax, returns ierr
  virtual int Apply(const TreeVector& x, TreeVector& b) const = 0;

  // Apply the inverse, x <-- A^-1 b, returns ierr
  virtual int ApplyInverse(const TreeVector& b, TreeVector& x) const = 0;
};

} // namespace Amanzi

#endif
