/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon
*/


//!  Timestep controller providing constant timestep size.

/*!

``TimestepControllerFixed`` is a simple timestep control mechanism which sets
a constant timestep size.  Note that the actual timestep size is given by the
minimum of PK's initial timestep sizes.

No parameters are required.

*/


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
      Errors::TimeStepCrash msg("Timestep failed: fixed time step size failed.");
      Exceptions::amanzi_throw(msg);
    }

    return dt;
  }

 protected:
  Teuchos::ParameterList plist_;
};

} // namespace Amanzi

#endif
