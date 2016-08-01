/*
  Time Integration

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon

  Simple timestep control based upon previous iteration count.
*/

#ifndef AMANZI_STANDARD_TIMESTEP_CONTROLLER_HH_
#define AMANZI_STANDARD_TIMESTEP_CONTROLLER_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "TimestepController.hh"

namespace Amanzi {

class TimestepControllerStandard : public TimestepController {

 public:
  TimestepControllerStandard(Teuchos::ParameterList& plist);

  // single method for timestep control
  double get_timestep(double dt, int iterations);

 protected:
  Teuchos::ParameterList plist_;

  int max_its_;
  int min_its_;
  double reduction_factor_;
  double increase_factor_;
  double max_dt_;
  double min_dt_;
};

} // namespace Amanzi

#endif
