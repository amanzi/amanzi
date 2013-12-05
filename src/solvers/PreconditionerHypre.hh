/*
  This is the Linear Solver component of the Amanzi code.

  License: BSD
  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)

  HYPRE preconditioner.
  Usage:
*/

#ifndef AMANZI_PRECONDITIONER_HYPRE_HH_
#define AMANZI_PRECONDITIONER_HYPRE_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_MultiVector.h"
#include "Epetra_RowMatrix.h"
#include "Ifpack.h"

#include "exceptions.hh"
#include "Preconditioner.hh"

namespace Amanzi {
namespace AmanziPreconditioners {

class PreconditionerHypre : public Preconditioner {
 public:
  PreconditionerHypre() {};
  ~PreconditionerHypre() {};

  void Init(const std::string& name, const Teuchos::ParameterList& list);
  void Update(const Teuchos::RCP<Epetra_RowMatrix>& A);
  void Destroy() {};

  int ApplyInverse(const Epetra_MultiVector& v, Epetra_MultiVector& hv);

 private:
  Teuchos::ParameterList list_;

  Teuchos::RCP<Ifpack_Hypre> IfpHypre_;
  double tol_, strong_threshold_;
  int nsmooth_, ncycles_, verbosity_;
};

}  // namespace AmanziPreconditioners
}  // namespace Amanzi



#endif
