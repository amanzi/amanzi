/*
  This is the Linear Solver component of the Amanzi code.

  License: BSD
  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)

  Incomplete LU preconditioner.
  Usage:
*/

#ifndef AMANZI_PRECONDITIONER_BLOCK_ILU_HH_
#define AMANZI_PRECONDITIONER_BLOCK_ILU_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_MultiVector.h"
#include "Epetra_RowMatrix.h"

#include "exceptions.hh"
#include "Preconditioner.hh"

namespace Amanzi {
namespace AmanziPreconditioners {

class PreconditionerBlockILU : public Preconditioner {
 public:
  PreconditionerBlockILU() {};
  ~PreconditionerBlockILU() {};

  void Init(const std::string& name, const Teuchos::ParameterList& list);
  void Update(const Teuchos::RCP<Epetra_RowMatrix>& A);
  void Destroy();

  void ApplyInverse(const Epetra_MultiVector& v, Epetra_MultiVector& hv);

 private:
  Teuchos::ParameterList list_;
  Teuchos::RCP<Ifpack_Preconditioner> IfpILU_;

  bool initialized_;
};

}  // namespace AmanziPreconditioners
}  // namespace Amanzi



#endif
