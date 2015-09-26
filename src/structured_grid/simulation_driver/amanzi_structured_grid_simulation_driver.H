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
 public:

  // v1 constructor
  explicit AmanziStructuredGridSimulationDriver(const std::string& xmlInFileName);

  // v2 constructor
  explicit AmanziStructuredGridSimulationDriver(xercesc::DOMDocument* input);

  // Destructor
  ~AmanziStructuredGridSimulationDriver();

  ReturnType Run(const MPI_Comm&               mpi_comm,
                 Amanzi::ObservationData&      output_observations);

  private:

  Teuchos::ParameterList* plist_;
};
#endif
