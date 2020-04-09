/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon (coonet@ornl.gov)
*/

//!

#ifndef AMANZI_TREEMATRIX_HH_
#define AMANZI_TREEMATRIX_HH_

#include "Teuchos_RCP.hpp"

namespace Amanzi {

class TreeVector;
class TreeVectorSpace;

class TreeMatrix {
 public:
  // Vector space of the Matrix's domain.
  virtual const TreeVectorSpace& DomaingetMap() const = 0;

  // Vector space of the Matrix's range.
  virtual const TreeVectorSpace& RangegetMap() const = 0;

  // Apply matrix, b <-- Ax, returns ierr
  virtual int apply(const TreeVector& x, TreeVector& b) const = 0;

  // Apply the inverse, x <-- A^-1 b, returns ierr
  virtual int applyInverse(const TreeVector& b, TreeVector& x) const = 0;
};

} // namespace Amanzi

#endif
