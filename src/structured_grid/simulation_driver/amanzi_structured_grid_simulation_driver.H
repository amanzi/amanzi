#ifndef _amanzi_structured_grid_simulation_driver_H
#define _amanzi_structured_grid_simulation_driver_H

#include "BoxLib.H"
#include "ParallelDescriptor.H"
#include "ObservationData.hh"
#include "Simulator.hh"

XERCES_CPP_NAMESPACE_USE

struct AmanziStructuredGridSimulationDriver : Amanzi::Simulator {
 public:

  // v2 constructor
  AmanziStructuredGridSimulationDriver(const std::string& input_file, DOMDocument* input,
                                       const std::string& ouput_prefix);

  // Destructor
  ~AmanziStructuredGridSimulationDriver();

  ReturnType Run(const Amanzi::Comm_ptr_type& comm,
                 Amanzi::ObservationData& output_observations);
};
#endif
