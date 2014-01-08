/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Fixed timestep size, with step size set by the PK's initial timestep size.
------------------------------------------------------------------------- */

#ifndef AMANZI_FIXED_TIMESTEP_CONTROLLER_HH_
#define AMANZI_FIXED_TIMESTEP_CONTROLLER_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "errors.hh"
#include "TimestepController.hh"

namespace Amanzi {

class TimestepControllerFixed : public TimestepController {

 public:
  TimestepControllerFixed(Teuchos::ParameterList& plist) : plist_(plist) {}

  // single method for timestep control
  double get_timestep(double dt, int iterations) {
    if (iterations < 0) {
      std::string msg = "Timestep failed: Time step crash";
      Errors::Message m(msg);
      Exceptions::amanzi_throw(m);
    }

    return dt;
  }

 protected:
  Teuchos::ParameterList plist_;
};

} // namespace

#endif
