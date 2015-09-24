#ifndef _amanzi_structured_grid_simulation_driver_H
#define _amanzi_structured_grid_simulation_driver_H

#include "BoxLib.H"
#include "ParallelDescriptor.H"
#include "Teuchos_ParameterList.hpp"
#include "ObservationData.hh"
#include "Simulator.hh"
#include "InputParser_Structured.H"

struct AmanziStructuredGridSimulationDriver
  : Amanzi::Simulator
{
  virtual ReturnType Run (const MPI_Comm&               mpi_comm,
                          Teuchos::ParameterList&       input_parameter_list,
                          Amanzi::ObservationData&      output_observations);
};
#endif
