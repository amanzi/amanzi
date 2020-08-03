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
namespace AmanziSolvers {

class PreconditionerDiagonal : public Preconditioner {
 public:
  virtual void InitializeInverse(Teuchos::ParameterList& plist) override final {};
  virtual void UpdateInverse() override final {}
  virtual void ComputeInverse() override final {
    AMANZI_ASSERT(h_.get());
    diagonal_ = Teuchos::rcp(new Epetra_Vector(h_->DomainMap()));
    h_->ExtractDiagonalCopy(*diagonal_);
  };

  virtual int ApplyInverse(const Epetra_Vector& v, Epetra_Vector& hv) const override final {
    AMANZI_ASSERT(diagonal_.get()); // Compute called
    return hv.ReciprocalMultiply(1.0, *diagonal_, v, 0.0); 
  }

  virtual int returned_code() const override final  { return 1; }
  virtual std::string returned_code_string() const override final { return "success"; }

 private:
  Teuchos::RCP<Epetra_Vector> diagonal_;
};

}  // namespace AmanziSolvers
}  // namespace Amanzi



#endif
