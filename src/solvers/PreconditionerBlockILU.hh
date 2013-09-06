/*
This is the Linear Solver component of the Amanzi code.
 
License: BSD
Authors: Konstantin Lipnikov (lipnikov@lanl.gov)

Incomplete LU preconditioner.
Usage: 
*/

#ifndef __PRECONDITIONER_BLOCK_ILU_HH__
#define __PRECONDITIONER_BLOCK_ILU_HH__

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_Vector.h"
#include "Epetra_FECrsMatrix.h"

#include "exceptions.hh"
#include "Preconditioner.hh"
 
namespace Amanzi {
namespace AmanziPreconditioners {

class PreconditionerBlockILU : public Preconditioner {
 public:
  PreconditionerBlockILU() {};
  ~PreconditionerBlockILU() {};

  void Init(const std::string& name, const Teuchos::ParameterList& list);
  void Update(Teuchos::RCP<Epetra_FECrsMatrix> A);
  void Destroy();

  void ApplyInverse(const Epetra_Vector& v, Epetra_Vector& hv);

 private:
  Teuchos::ParameterList list_;
  Teuchos::RCP<Ifpack_Preconditioner> IfpILU_;

  bool initialized;
};

}  // namespace AmanziPreconditioners
}  // namespace Amanzi



#endif


