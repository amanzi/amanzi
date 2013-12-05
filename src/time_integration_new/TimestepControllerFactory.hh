/*
  Factory for timestep control.

  License: BSD
  Authors: Ethan Coon (ecoon@lanl.gov)

*/

#ifndef AMANZI_TS_CONTROLLER_FACTORY_HH_
#define AMANZI_TS_CONTROLLER_FACTORY_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "errors.hh"
#include "TimestepController.hh"

namespace Amanzi {

struct TimestepControllerFactory {
 public:
  Teuchos::RCP<TimestepController>
  Create(const Teuchos::ParameterList& prec_list);
};

}  // namespace Amanzi

#endif
