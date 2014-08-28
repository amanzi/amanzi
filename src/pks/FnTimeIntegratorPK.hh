/*
  License: see $AMANZI_DIR/COPYRIGHT
  Authors: Konstantin Lipnikov, Ethan Coon

  This is a purely virtual base class for process kernels which use
  time integrators.

*/

#ifndef AMANZI_FN_TIME_INTEGRATOR_PK_HH_
#define AMANZI_FN_TIME_INTEGRATOR_PK_HH_

#include "Teuchos_RCP.hpp"

#include "BDFFnBase.hh"
#include "PK.hh"

namespace Amanzi {

class TreeVector;

class FnTimeIntegratorPK : public PK, public Amanzi::BDFFnBase<TreeVector> {};

}  // namespace Amanzi

#endif
