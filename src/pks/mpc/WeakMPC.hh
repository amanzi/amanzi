#ifndef _MPC_HPP_
#define _MPC_HPP_

#include <vector>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Epetra_MpiComm.h"

#include "State.hh"
#include "Vector.hh"
#include "TreeVector.hh"

#include "PK.hh"
#include "PK_Factory.hh"
#include "MPC.hh"

namespace Amanzi {

class WeakMPC : public Teuchos::VerboseObject<WeakMPC>, public MPC {

public:
  WeakMPC(Teuchos::ParameterList& mpc_plist,
          Teuchos::RCP<State>& S, Teuchos::RCP<TreeVector>& soln);

  ~WeakMPC() {};

  bool advance(double dt, Teuchos::RCP<const State> &S0,
               Teuchos::RCP<State> &S1, Teuchos::RCP<TreeVector> &solution);
};
} // close namespace Amanzi

#endif
