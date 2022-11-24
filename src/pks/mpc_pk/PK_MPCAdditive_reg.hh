/*
  This is the mpc_pk component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon

  Implementation for the derived PK_MPCStrong class.  Is both a PK and a Model
  Evalulator, providing needed methods for BDF time integration of the coupled
  system.

  Completely automated and generic to any sub PKs, this uses a block diagonal
  preconditioner.

  See additional documentation in the base class src/pks/mpc_pk/PK_MPC.hh
*/

#include "PK_MPCAdditive.hh"

namespace Amanzi {

template <>
RegisteredPKFactory<PK_MPCStrong<PK_Implicit>> PK_MPCStrong<PK_Implicit>::reg_("additive MPC");

} // namespace Amanzi
