/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Jeffrey Johnson <jnjohnson@lbl.gov>
*/

/*
  Simulator

*/

#ifndef AMANZI_SIMULATOR_FACTORY_HH_
#define AMANZI_SIMULATOR_FACTORY_HH_

#include "Simulator.hh"

// This factory method creates a Simulator object that can run a simulation
// for the input contained by the file with the given name.
namespace Amanzi {
namespace SimulatorFactory {

std::unique_ptr<Simulator> Create(const std::string& input_filename,
                                  const std::string& output_prefix);

} // namespace SimulatorFactory
} // namespace Amanzi

#endif
