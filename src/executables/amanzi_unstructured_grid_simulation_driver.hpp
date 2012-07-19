
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "ObservationData.H"
#include "Simulator.H"


// for printing the native parameter list to xml we need to 
// make sure that double gets printed with a sufficient number
// of digits
namespace Teuchos {
  template<>
  class ToStringTraits<double> {
  public:
    static std::string toString (const double& t) {
      std::ostringstream os;
      os.setf (std::ios::scientific);
      os.precision (17);
      os << t;
      return os.str();
    }
  };
}


struct AmanziUnstructuredGridSimulationDriver
  : Amanzi::Simulator, public Teuchos::VerboseObject<AmanziUnstructuredGridSimulationDriver>
{
  virtual ReturnType Run (const MPI_Comm&               mpi_comm,
                          Teuchos::ParameterList&       input_parameter_list,
                          Amanzi::ObservationData&      output_observations);
};
