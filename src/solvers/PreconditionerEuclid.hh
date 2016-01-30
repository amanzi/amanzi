/*
  Solvers

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
           Ethan Coon (ecoon@lanl.gov)

  HYPRE Euclid parallel ILU preconditioner.
*/

#ifndef AMANZI_PRECONDITIONER_EUCLID_HH_
#define AMANZI_PRECONDITIONER_EUCLID_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_MultiVector.h"
#include "Epetra_RowMatrix.h"
#include "Ifpack.h"

#include "exceptions.hh"
#include "Preconditioner.hh"

namespace Amanzi {
namespace AmanziPreconditioners {

class PreconditionerEuclid : public Preconditioner {
 public:
  PreconditionerEuclid() {};
  ~PreconditionerEuclid() {};

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
