#include "mpi.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "ObservationData.hh"
#include "Simulator.hh"


struct AmanziUnstructuredGridSimulationDriver
    : Amanzi::Simulator, 
      public Teuchos::VerboseObject<AmanziUnstructuredGridSimulationDriver>
{
 public:
  // v1 constructor
  // explicit AmanziUnstructuredGridSimulationDriver(const std::string& xmlInFileName);

  // v2 constructor
  AmanziUnstructuredGridSimulationDriver(const std::string& xmlInFileName,
                                         xercesc::DOMDocument* input);

  // Destructor.
  ~AmanziUnstructuredGridSimulationDriver();

  ReturnType Run(const MPI_Comm& mpi_comm,
                 Amanzi::ObservationData& observations_data);

 private:
  // Read our parameter list.
  void ReadParameterList();
  Teuchos::ParameterList* plist_;
};
