/*
This is the Linear Solver component of the Amanzi code.
 
License: BSD
Authors: Konstantin Lipnikov (lipnikov@lanl.gov)

Conjugate gradient method.
Usage: 
*/

#ifndef __PRECONDITIONER_ML_HH__
#define __PRECONDITIONER_ML_HH__

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_Vector.h"
#include "Epetra_FECrsMatrix.h"
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
  void Update(Teuchos::RCP<Epetra_FECrsMatrix> A);
  void Destroy();

  void ApplyInverse(const Epetra_Vector& v, Epetra_Vector& hv);

 private:
  Teuchos::ParameterList list_;
  Teuchos::RCP<ML_Epetra::MultiLevelPreconditioner> ML_;

  bool initialized;
};

}  // namespace AmanziPreconditioners
}  // namespace Amanzi



#endif


