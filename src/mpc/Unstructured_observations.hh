#ifndef _UNSTRUCTURED_OBSERVATIONS_HH_
#define _UNSTRUCTURED_OBSERVATIONS_HH_

#include <map>

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Array.hpp"

#include "ObservationData.hh"
#include "State.hh"
#include "io_event.hh"

#include "TimeStepManager.hh"
#include "VerboseObject.hh"


namespace Amanzi {

class Unstructured_observations {
 public:
  struct Observable : public IOEvent {
    Observable(std::string variable_,
               std::string region_,
               std::string functional_,
               Teuchos::ParameterList& plist,
               Epetra_MpiComm* comm):
        variable(variable_), region(region_),
        functional(functional_), plist_(plist),
        IOEvent(plist, comm)
    {
      ReadParameters_();
    }

    std::string variable;
    std::string region;
    std::string functional;
    const Teuchos::ParameterList& plist_;
  };

  Unstructured_observations(Teuchos::ParameterList obs_list,
                            Amanzi::ObservationData& observation_data,
			    Epetra_MpiComm* comm);
  
  void RegisterComponentNames(std::vector<std::string> comp_names) {
    comp_names_ = comp_names;
  }

  void MakeObservations(State& state);

  bool DumpRequested(const int);
  bool DumpRequested(const double);
  bool DumpRequested(const int, const double);
     
  void RegisterWithTimeStepManager(const Teuchos::Ptr<TimeStepManager>& tsm);

  void FlushObservations();

 protected:
  VerboseObject* vo_;
  
 private:
  int rank_;
  Teuchos::ParameterList obs_list_;
  Amanzi::ObservationData& observation_data_;
  std::map<std::string, Observable> observations;
  std::vector<std::string> comp_names_;
  std::map<std::string, double> drawdown_;
};

}  // namespace Amanzi

#endif
