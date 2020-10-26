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

This is provided when using the `"preconditioning method`"=`"diagonal`" in the
`Preconditioner`_ spec.

No parameters are required.

*/


#ifndef AMANZI_PRECONDITIONER_DIAGONAL_HH_
#define AMANZI_PRECONDITIONER_DIAGONAL_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
//#include "Epetra_MultiVector.h"
//#include "Epetra_RowMatrix.h"

#include "exceptions.hh"
#include "Preconditioner.hh"

namespace Amanzi {
namespace AmanziSolvers {

class PreconditionerDiagonal : public Preconditioner {
 public:
  virtual void set_inverse_parameters(Teuchos::ParameterList& plist) override final {};
  virtual void initializeInverse() override final {}
  virtual void computeInverse() override final {
    diagonal_ = Teuchos::rcp(new Vector_type(h_->getRowMap()));
    diagonal_->putScalar(0.);
    h_->getLocalDiagCopy(*diagonal_);
    diagonal_->reciprocal(*diagonal_);
  };

  virtual int applyInverse(const Vector_type& v, Vector_type& hv) const override final {
    AMANZI_ASSERT(diagonal_.get()); // Compute called
    hv.elementWiseMultiply(1., v, *diagonal_, 0.);
    return 0;
  }

  virtual int returned_code() const override final  { return returned_code_; }
  virtual std::string returned_code_string() const override final {
    if (returned_code_ == 0) return "success";
    return "failed ReciprocalMultiply()";
  }

  virtual void update(const Teuchos::RCP<Matrix_type>&) override final {}; 


 private:
  Teuchos::RCP<Vector_type> diagonal_;
  mutable int returned_code_;
};

}  // namespace AmanziSolvers
}  // namespace Amanzi



#endif
