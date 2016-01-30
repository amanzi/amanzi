/*
  Solvers

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  ML preconditioner.
*/

#ifndef AMANZI_PRECONDITIONER_ML_HH_
#define AMANZI_PRECONDITIONER_ML_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_MultiVector.h"
#include "Epetra_RowMatrix.h"
#include "ml_MultiLevelPreconditioner.h"

#include "exceptions.hh"
#include "Preconditioner.hh"

namespace Amanzi {
namespace AmanziPreconditioners {

class PreconditionerML : public Preconditioner {
 public:
  PreconditionerML() {};
  ~PreconditionerML() {};

  void Init(const std::string& name, const Teuchos::ParameterList& list);
  void Update(const Teuchos::RCP<Epetra_RowMatrix>& A);
  void Destroy();

  int ApplyInverse(const Epetra_MultiVector& v, Epetra_MultiVector& hv);

  int returned_code() { return returned_code_; }

 private:
  Teuchos::ParameterList list_;
  Teuchos::RCP<ML_Epetra::MultiLevelPreconditioner> ML_;

  bool initialized_;
  int returned_code_;
};

}  // namespace AmanziPreconditioners
}  // namespace Amanzi



#endif
