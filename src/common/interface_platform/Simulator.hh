#ifndef SIMULATOR_H
#define SIMULATOR_H

#include "map"
#include "vector"

#include "ObservationData.hh"
#include "Teuchos_ParameterList.hpp"

namespace Amanzi
{
  struct Simulator
  {
    enum ReturnType {SUCCESS, FAIL, NUM_RETURN_TYPES};
    
    virtual Amanzi::Simulator::ReturnType Run(const MPI_Comm&               mpi_comm,
                                              Teuchos::ParameterList&       input_parameter_list,
                                              Amanzi::ObservationData&      output_observations) = 0;
  };
} // end namespace amanzi

#endif
