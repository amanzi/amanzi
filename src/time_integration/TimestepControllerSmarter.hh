/*
  Time Integration

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon

  Slightly smarter timestep control based upon a history of previous timesteps.
*/

#ifndef AMANZI_SMARTER_TIMESTEP_CONTROLLER_HH_
#define AMANZI_SMARTER_TIMESTEP_CONTROLLER_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "TimestepController.hh"

namespace Amanzi {

class TimestepControllerSmarter : public TimestepController {

 public:
  TimestepControllerSmarter(Teuchos::ParameterList& plist);

  // single method for timestep control
  double get_timestep(double dt, int iterations);

 protected:
  Teuchos::ParameterList plist_;

  int max_its_;
  int min_its_;

  double reduction_factor_;

  double increase_factor_;
  double increase_factor0_;
  double max_increase_factor_;
  int count_increased_before_increase_;
  int successive_increases_;

  double max_dt_;
  double min_dt_;

  int last_fail_;
  int growth_wait_after_fail_;
  int growth_wait_after_fail0_;
};

} // namespace

#endif
