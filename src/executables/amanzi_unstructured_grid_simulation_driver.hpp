
#include "Teuchos_ParameterList.hpp"
#include "ObservationData.H"
#include "Simulator.H"

struct AmanziUnstructuredGridSimulationDriver
  : Simulator
{
  virtual ReturnType Run (const MPI_Comm&               mpi_comm,
                          const Teuchos::ParameterList& input_parameter_list,
                          ObservationData&              output_observations);
};
