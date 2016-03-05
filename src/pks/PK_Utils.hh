/*
  This is the process kernel component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Miscalleneous collection of simple non-member functions.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_PK_UTILS_HH_
#define AMANZI_PK_UTILS_HH_

#include "State.hh"

namespace Amanzi {

// Averages permeability tensor in horizontal direction.
void PKUtils_CalculatePermeabilityFactorInWell(
    const Teuchos::Ptr<State>& S, Teuchos::RCP<Epetra_Vector>& Kxy);

}  // namespace Amanzi

#endif
