/*
  This is the mpc_pk component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov

  Process kernel for coupling Flow PK with Energy PK.
*/

#include "FlowEnergy_PK.hh"
#include "MPCStrong.hh"

namespace Amanzi {

// -----------------------------------------------------------------------------
// Constructor
// -----------------------------------------------------------------------------
FlowEnergy_PK::FlowEnergy_PK(Teuchos::ParameterList& pk_tree,
    const Teuchos::RCP<Teuchos::ParameterList>& global_list,
    const Teuchos::RCP<State>& S,
    const Teuchos::RCP<TreeVector>& soln) :
    Amanzi::MPCStrong<FnTimeIntegratorPK>(pk_tree, global_list, S, soln) { 
}

}  // namespace Amanzi

