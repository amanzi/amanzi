/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------

ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon (ecoon@lanl.gov)

Interface for a Matrix that acts on Epetra_Vector.

This is the local Amanzi version, not a true Epetra product.  This
takes an Epetra_RowMatrix and wraps it with the needed methods to be
used with solvers.
------------------------------------------------------------------------- */

#ifndef AMANZI_EPETRAMATRIX_HH_
#define AMANZI_EPETRAMATRIX_HH_

#include "Teuchos_RCP.hpp"

namespace Amanzi {

class EpetraMatrix {
 public:
  // Vector space of the Matrix's domain.
  virtual const Epetra_BlockMap& DomainMap() const = 0;

  // Vector space of the Matrix's range.
  virtual const Epetra_BlockMap& RangeMap() const = 0;

  // Apply matrix, b <-- Ax, returns ierr
  virtual int Apply(const Epetra_Vector& x, Epetra_Vector& b) const = 0;

  // Apply the inverse, x <-- A^-1 b, returns ierr
  virtual int ApplyInverse(const Epetra_Vector& b, Epetra_Vector& x) const = 0;
};

} // namespace Amanzi

#endif
