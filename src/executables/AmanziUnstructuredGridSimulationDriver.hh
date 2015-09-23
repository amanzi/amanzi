#include "mpi.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "ObservationData.hh"
#include "Simulator.hh"


struct AmanziUnstructuredGridSimulationDriver
    : Amanzi::Simulator, 
      public Teuchos::VerboseObject<AmanziUnstructuredGridSimulationDriver>
{
  explicit AmanziUnstructuredGridSimulationDriver(const Teuchos::ParameterList& input_parameter_list);

  virtual ReturnType Run(const MPI_Comm& mpi_comm,
                         Amanzi::ObservationData& output_observations);

  private:

  Teuchos::ParameterList param_list_;
};
