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

This is provided when using the `"preconditioner type`"=`"identity`" in the
`Preconditioner`_ spec.

No parameters are required.

*/


#ifndef AMANZI_PRECONDITIONER_IDENTITY_HH_
#define AMANZI_PRECONDITIONER_IDENTITY_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_MultiVector.h"
#include "Epetra_RowMatrix.h"

#include "exceptions.hh"
#include "Preconditioner.hh"

namespace Amanzi {
namespace AmanziPreconditioners {

class PreconditionerIdentity : public Preconditioner {
 public:
  PreconditionerIdentity() {};
  ~PreconditionerIdentity() {};

  void Init(const std::string& name, const Teuchos::ParameterList& list) {};
  void Update(const Teuchos::RCP<Epetra_RowMatrix>& A) {};
  void Destroy() {};

  int ApplyInverse(const Epetra_MultiVector& v, Epetra_MultiVector& hv) {
    hv = v;
    return 0;
  }

  int returned_code() { return 0; }
};

}  // namespace AmanziPreconditioners
}  // namespace Amanzi



#endif
