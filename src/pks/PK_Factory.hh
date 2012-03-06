/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Interface for the process kernel factory.
------------------------------------------------------------------------- */

#ifndef PKS_PK_FACTORY_HH_
#define PKS_PK_FACTORY_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "State.hh"
#include "PK.hh"

namespace Amanzi {

class PK_Factory {

public:

  Teuchos::RCP<PK> create_pk(Teuchos::ParameterList plist, Teuchos::RCP<State>& S,
                             Teuchos::RCP<TreeVector>& soln);
};
}

#endif
