#ifndef __DARCY_MATVEC_H__
#define __DARCY_MATVEC_H__

#include "Teuchos_RCP.hpp"
#include "Epetra_Operator.h"

#include "DarcyProblem.hpp"

namespace Amanzi {
namespace AmanziFlow {

class DarcyMatvec : public Epetra_Operator {
 public:
  DarcyMatvec(DarcyProblem *problem);
  ~DarcyMatvec() {}

  int Apply(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const;

  // The remaining methods return errors or hardwired results.
  // It's not expected that any of them should be called.

  int ApplyInverse(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const { return 1; }

  double NormInf() const { return 0.0; }

  const char* Label() const { return label_; }

  bool UseTranspose() const { return false; }

  int SetUseTranspose(bool) { return 1; }

  bool HasNormInf() const { return false; }

  const Epetra_Comm& Comm() const { return problem_->Comm(); }

  const Epetra_Map& OperatorDomainMap() const { return problem_->Map(); }
  const Epetra_Map& OperatorRangeMap() const { return problem_->Map(); }

 private:
  char* label_;

  DarcyProblem *problem_;
};

}  // namespace AmanziFlow
}  // namespace Amanzi

#endif
