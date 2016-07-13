/*
  Simulator 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Jeffrey Johnson <jnjohnson@lbl.gov>
*/

#ifndef AMANZI_SIMULATOR_FACTORY_H
#define AMANZI_SIMULATOR_FACTORY_H

#include "Simulator.hh"

// This factory method creates a Simulator object that can run a simulation 
// for the input contained by the file with the given name.
namespace Amanzi {
namespace SimulatorFactory {

Simulator* Create(const std::string& input_filename);

} // namespace SimulatorFactory
} // namespace Amanzi

#endif
