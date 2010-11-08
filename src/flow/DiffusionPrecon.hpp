#ifndef __DIFFUSIONPRECON_H__
#define __DIFFUSIONPRECON_H__

#include "Teuchos_RCP.hpp"
#include "Epetra_Operator.h"

#include "DiffusionMatrix.hpp"

#include "ml_MultiLevelPreconditioner.h"

class DiffusionPrecon : public Epetra_Operator
{
public:
  DiffusionPrecon(Teuchos::RCP<DiffusionMatrix> &matrix, const Epetra_Map &map);
  ~DiffusionPrecon() { };

  void Compute();

  int ApplyInverse(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const;

  // The remaining methods are required by the Epetra_Operator interface but are not
  // expected to be used, and are defined to return errors or hardwired results.

  int Apply(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const { return 1; }

  double NormInf() const { return 0.0; }

  const char* Label() const { return label; }

  bool UseTranspose() const { return false; }

  int SetUseTranspose(bool) { return 1; }

  bool HasNormInf() const { return false; }

  const Epetra_Comm& Comm() const { return D->Comm(); }

  const Epetra_Map& OperatorDomainMap() const { return map_; }

  const Epetra_Map& OperatorRangeMap() const { return map_; }

private:
  // The diffusion matrix
  Teuchos::RCP<DiffusionMatrix> D;

  Epetra_Map map_;

  // ML preconditioner for the Schur complement system.
  ML_Epetra::MultiLevelPreconditioner *MLprec;

  char *label;
};

#endif
