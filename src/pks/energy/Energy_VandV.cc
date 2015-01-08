/*
  This is the energy component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <set>

#include "errors.hh"
#include "OperatorDefs.hh"

#include "Energy_PK.hh"

namespace Amanzi {
namespace Energy {

/* ****************************************************************
* Construct default state for unit tests. It completes
* initialization of all missed objects in the state.
**************************************************************** */
void Energy_PK::InitializeFields()
{
  // set popular default values
  /*
  if (!S_->GetField("fluid_density", passwd_)->initialized()) {
    *(S_->GetScalarData("fluid_density", passwd_)) = 1000.0;
    S_->GetField("fluid_density", passwd_)->set_initialized();
  }
  */

  if (!S_->GetField("temperature", passwd_)->initialized()) {
    S_->GetFieldData("temperature", passwd_)->PutScalar(300.0);
    S_->GetField("temperature", passwd_)->set_initialized();
  }
}

}  // namespace Energy
}  // namespace Amanzi

