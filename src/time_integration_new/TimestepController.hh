/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Virtual interface for timestep adaptation.
------------------------------------------------------------------------- */

#ifndef AMANZI_TIMESTEP_CONTROLLER_HH_
#define AMANZI_TIMESTEP_CONTROLLER_HH_

#include "Teuchos_RCP.hpp"

namespace Amanzi {

class TimestepController {

public:
  // virtual destructor
  virtual ~TimestepController() {};

  // single method for timestep control
  virtual double get_timestep(double dt, int iterations) = 0;

};

} // namespace

#endif
