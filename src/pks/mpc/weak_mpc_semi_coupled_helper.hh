#ifndef WEAK_MPC_SEMI_COUPLED_HELPER_HH_
#define WEAK_MPC_SEMI_COUPLED_HELPER_HH_

#include "State.hh"

namespace Amanzi{

void
UpdateIntermediateStateParameters(Teuchos::RCP<Amanzi::State>& S_next_, Teuchos::RCP<Amanzi::State>& S_inter_, int id);

void
UpdateNextStateParameters(Teuchos::RCP<Amanzi::State>& S_next_, Teuchos::RCP<Amanzi::State>& S_inter_, int id);



}

#endif
