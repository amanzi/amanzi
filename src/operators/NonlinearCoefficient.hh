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

namespace Amanzi {
namespace Operators {

class NonlinearCoefficient {
 public:
  NonlinearCoefficient(const Teuchos::RCP<Epetra_Vector> cdata) 
      : cdata_(cdata) {};
  ~NonlinearCoefficient() {};

  // main methods
  void UpdateCellValues() {};
  void UpwindCellValues();

  // access
  const Teuchos::RCP<Epetra_Vector> cdata() { return cdata_; }

 private:
  const Teuchos::RCP<Epetra_Vector> cdata_;
  Teuchos::RCP<Epetra_Vector> fdata_;
};

}  // namespace Operators
}  // namespace Amanzi

#endif

