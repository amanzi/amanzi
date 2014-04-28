/*
  This is the Linear Solver component of the Amanzi code.

  License: BSD
  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)

  HYPRE preconditioner.
  Usage:
*/

#ifndef AMANZI_PRECONDITIONER_BOOMERAMG_HH_
#define AMANZI_PRECONDITIONER_BOOMERAMG_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_MultiVector.h"
#include "Epetra_RowMatrix.h"
#include "Ifpack.h"

#include "exceptions.hh"
#include "Preconditioner.hh"

namespace Amanzi {
namespace AmanziPreconditioners {

class PreconditionerBoomerAMG : public Preconditioner {
 public:
  PreconditionerBoomerAMG() {};
  ~PreconditionerBoomerAMG() {};

  void Init(const std::string& name, const Teuchos::ParameterList& list);
  void Update(const Teuchos::RCP<Epetra_RowMatrix>& A);
  void Destroy() {};

  int ApplyInverse(const Epetra_MultiVector& v, Epetra_MultiVector& hv);

  int returned_code() { return returned_code_; }

 private:
  Teuchos::ParameterList plist_;
  std::vector<Teuchos::RCP<FunctionParameter> > funcs_;
  
  int returned_code_;
  Teuchos::RCP<Ifpack_Hypre> IfpHypre_;
};

}  // namespace AmanziPreconditioners
}  // namespace Amanzi



#endif
