/*
  Operators 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_OPERATOR_AUDIT_HH_
#define AMANZI_OPERATOR_AUDIT_HH_

#include "Epetra_CrsMatrix.h"

#include "Op.hh"


namespace Amanzi {
namespace Operators {

int CheckMatrixSymmetry(Teuchos::RCP<Epetra_CrsMatrix> A);
int CheckMatrixCoercivity(Teuchos::RCP<Epetra_CrsMatrix> A);

class OperatorAudit {
 public:
  OperatorAudit(Teuchos::RCP<const Op> op) : op_(op) {};
  ~OperatorAudit() {};

  // main members
  int CheckSpectralBounds(int schema);

 private:
  void OrderByIncrease_(int n, double* mem);

 private:
  Teuchos::RCP<const Op> op_;
};

}  // namespace Operators
}  // namespace Amanzi

#endif
