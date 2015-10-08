#ifndef SIMULATORFACTORY_H
#define SIMULATORFACTORY_H

#include "xercesc/util/PlatformUtils.hpp"
#include "Simulator.hh"

// This factory method creates a Simulator object that can run a simulation 
// for the input contained by the file with the given name.
namespace Amanzi
{

  namespace SimulatorFactory
  {
    Simulator* Create(const std::string& input_filename);
  }

} // end namespace amanzi

#endif
