#ifndef SIMULATOR_H
#define SIMULATOR_H

#include "map"
#include "vector"

#include "ObservationData.H"
#include "Teuchos_ParameterList.hpp"

struct Simulator
{
  enum ReturnType {SUCCESS, FAIL, NUM_RETURN_TYPES};

  virtual ReturnType Run (const MPI_Comm&               mpi_comm,
                          const Teuchos::ParameterList& input_parameter_list,
                          ObservationData&              output_observations) = 0;
};

#endif
