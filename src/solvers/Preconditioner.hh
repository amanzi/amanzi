/*
  Solvers

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)

  Base class for preconditioners.
*/

#ifndef AMANZI_PRECONDITIONER_HH_
#define AMANZI_PRECONDITIONER_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_MultiVector.h"
#include "Epetra_RowMatrix.h"

#include "exceptions.hh"

namespace Amanzi {
namespace AmanziPreconditioners {

class Preconditioner {
 public:
  virtual ~Preconditioner() = default;

  // Initializes the solver with provided parameters.
  // This need not be called by preconditioners created using the factory.
  virtual void Init(const std::string& name,
                    const Teuchos::ParameterList& list) = 0;

  // Rebuild the preconditioner using the given matrix A.
  virtual void Update(const Teuchos::RCP<Epetra_RowMatrix>& A) = 0;

  // Destroy the preconditioner and auxiliary data structures.
  virtual void Destroy() = 0;

  // Apply the preconditioner.
  virtual int ApplyInverse(const Epetra_MultiVector& v,
                           Epetra_MultiVector& hv) = 0;

  virtual int returned_code() = 0;
};

}  // namespace AmanziPreconditioners
}  // namespace Amanzi


#endif


