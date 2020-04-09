/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Konstantin Lipnikov (lipnikov@lanl.gov)
      Ethan Coon (coonet@ornl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>

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

class PreconditionerEuclid
  : public Preconditioner<Epetra_RowMatrix, Epetra_MultiVector> {
 public:
  PreconditionerEuclid(){};
  ~PreconditionerEuclid(){};

  void
  Init(const std::string& name, const Teuchos::ParameterList& list) override;
  void Update(const Teuchos::RCP<const Epetra_RowMatrix>& A) override;
  void Destroy() override{};

  int ApplyInverse(const Epetra_MultiVector& v,
                   Epetra_MultiVector& hv) const override;

  int returned_code() override { return returned_code_; }

 private:
  Teuchos::ParameterList plist_;
  std::vector<Teuchos::RCP<FunctionParameter>> funcs_;

  mutable int returned_code_;
  Teuchos::RCP<Ifpack_Hypre> IfpHypre_;
};

} // namespace AmanziPreconditioners
} // namespace Amanzi


#endif
