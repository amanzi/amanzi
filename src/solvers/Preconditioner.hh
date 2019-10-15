/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Konstantin Lipnikov (lipnikov@lanl.gov)
*/

//! Base class for preconditioners.

/*!

``[preconditioner-typed-spec]``

* `"preconditioner type`" ``[string]`` Iterative method to be used.
* `"_preconditioner_type_ parameters`" ``[_preconditioner_type_-spec]``
  Parameters associated with the requested preconditioner.

 */

#ifndef AMANZI_PRECONDITIONER_HH_
#define AMANZI_PRECONDITIONER_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_MultiVector.h"
#include "Epetra_RowMatrix.h"

#include "exceptions.hh"

namespace Amanzi {
namespace AmanziPreconditioners {

template <class Matrix, class Vector>
class Preconditioner {
 public:
  virtual ~Preconditioner() = default;

  // Initializes the solver with provided parameters.
  // This need not be called by preconditioners created using the factory.
  virtual void
  Init(const std::string& name, const Teuchos::ParameterList& list) = 0;

  // Rebuild the preconditioner using the given matrix A.
  virtual void Update(const Teuchos::RCP<const Matrix>& A) = 0;

  // Destroy the preconditioner and auxiliary data structures.
  virtual void Destroy() = 0;

  // Apply the preconditioner.
  virtual int applyInverse(const Vector& v, Vector& hv) const = 0;

  virtual int returned_code() = 0;
};

} // namespace AmanziPreconditioners
} // namespace Amanzi


#endif
