/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/
//! Identity as a preconditioner.

/*!

Simply copies the input vector to the output -- uses the Identity matrix as a
preconditioner.

This is provided when using the `"preconditioning method`"=`"identity`" in the
`Preconditioner`_ spec.

No parameters are required.

*/


#ifndef AMANZI_PRECONDITIONER_IDENTITY_HH_
#define AMANZI_PRECONDITIONER_IDENTITY_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
//#include "Epetra_MultiVector.h"
//#include "Epetra_RowMatrix.h"

#include "exceptions.hh"
#include "Inverse.hh"

namespace Amanzi {
namespace AmanziSolvers {

template<class Matrix,
         class Preconditioner=Matrix,
         class Vector=typename Matrix::Vector_t,
         class VectorSpace=typename Vector::VectorSpace_t>
class PreconditionerIdentity :
      public Inverse<Matrix,Preconditioner,Vector,VectorSpace> {
 public:
  PreconditionerIdentity() {};
  ~PreconditionerIdentity() {};

  virtual void set_inverse_parameters(Teuchos::ParameterList& list) override final {};
  virtual void initializeInverse() override final {};
  virtual void computeInverse() override final {};

  virtual int applyInverse(const Vector& v, Vector& hv) const override final {
    hv.assign(v);
    return 0;
  }

  virtual int returned_code() const override final { return 0; }
  virtual std::string returned_code_string() const override {
    return "PreconditionerIdentity was applied.";
  }
};

}  // namespace AmanziSolvers
}  // namespace Amanzi



#endif
