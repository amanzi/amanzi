/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/
//! Diagonal preconditioner.

/*!

Simply applys the pointwise inverse of the diagonal of the matrix as an
extremely cheap matrix.

This is provided when using the `"preconditioner type`"=`"diagonal`" in the
`Preconditioner`_ spec.

No parameters are required.

*/


#ifndef AMANZI_PRECONDITIONER_DIAGONAL_HH_
#define AMANZI_PRECONDITIONER_DIAGONAL_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_MultiVector.h"
#include "Epetra_RowMatrix.h"

#include "exceptions.hh"
#include "Preconditioner.hh"

namespace Amanzi {
namespace AmanziPreconditioners {

class PreconditionerDiagonal : public Preconditioner {
 public:
  PreconditionerDiagonal() {};
  ~PreconditionerDiagonal() {};

  void Init(const std::string& name, const Teuchos::ParameterList& plist) {};
  void Update(const Teuchos::RCP<Epetra_RowMatrix>& A) {
    diagonal_ = Teuchos::rcp(new Epetra_Vector(A->Map()));
    A->ExtractDiagonalCopy(*diagonal_);
  };
  void Destroy() {};

  int ApplyInverse(const Epetra_MultiVector& v, Epetra_MultiVector& hv) {
    return hv.ReciprocalMultiply(1.0, *diagonal_, v, 0.0); 
  }

  int returned_code() { return 0; }

 private:
  Teuchos::RCP<Epetra_Vector> diagonal_;
};

}  // namespace AmanziPreconditioners
}  // namespace Amanzi



#endif
