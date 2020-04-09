/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/

//! Timestep controller which loads a timestep history from file.

/*!

``TimestepControllerFromFile`` loads a timestep history from a file, then
advances the step size with those values.  This is mostly used for testing
purposes, where we need to force the same timestep history as previous runs to
do regression testing.  Otherwise even machine roundoff can eventually alter
number of iterations enough to alter the timestep history, resulting in
solutions which are enough different to cause doubt over their correctness.

* `"file name`" ``[string]`` Path to hdf5 file containing timestep information.
* `"timestep header`" ``[string]`` Name of the dataset containing the history of
timestep sizes.

*/


#ifndef AMANZI_FROMFILE_TIMESTEP_CONTROLLER_HH_
#define AMANZI_FROMFILE_TIMESTEP_CONTROLLER_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "errors.hh"
#include "TimestepController.hh"

namespace Amanzi {

class TimestepControllerFromFile : public TimestepController {
 public:
  TimestepControllerFromFile(Teuchos::ParameterList& plist);

  // single method for timestep control
  double get_timestep(double dt, int iterations);

 protected:
  std::vector<double> dt_history_;
  int current_;
};

} // namespace Amanzi

#endif
