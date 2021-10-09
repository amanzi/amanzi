/*
  MPC PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov

  Sequential coupling of shallow water and solute transport.
*/

#include "ShallowWaterTransport_PK.hh"

namespace Amanzi {

/* ******************************************************************
* Setup of PK
****************************************************************** */
void ShallowWaterTransport_PK::Setup(const Teuchos::Ptr<State>& S)
{
  PK_MPCWeak::Setup(S);
}

}  // namespace Amanzi

