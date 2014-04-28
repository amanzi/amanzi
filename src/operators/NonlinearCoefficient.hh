/*
  This is the operators component of the Amanzi code. 

  Copyright 2010-2013 held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_NONLINEAR_COEFFICIENT_HH_
#define AMANZI_NONLINEAR_COEFFICIENT_HH_

#include <string>
#include <vector>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "CompositeVector.hh"

namespace Amanzi {
namespace Operators {

class NonlinearCoefficient {
 public:
  NonlinearCoefficient() {}; 
  ~NonlinearCoefficient() {};

  // main methods
  virtual void UpdateValues(const CompositeVector& u) {};
  virtual void UpdateDerivatives(const CompositeVector& u) {};

  // access
  const Teuchos::RCP<Epetra_Vector> cvalues() { return cvalues_; }
  const Teuchos::RCP<Epetra_Vector> fvalues() { return fvalues_; }

  const Teuchos::RCP<Epetra_Vector> fderivatives() { return fderivatives_; }

 public:
  Teuchos::RCP<Epetra_Vector> cvalues_;
  Teuchos::RCP<Epetra_Vector> fvalues_;

  Teuchos::RCP<Epetra_Vector> fderivatives_;
};

}  // namespace Operators
}  // namespace Amanzi

#endif

