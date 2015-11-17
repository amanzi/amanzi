#ifndef _UNSTRUCTURED_OBSERVATIONS_HH_
#define _UNSTRUCTURED_OBSERVATIONS_HH_

#include <map>

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Array.hpp"

#include "ObservationData.hh"
#include "State.hh"
#include "IOEvent.hh"

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
        IOEvent(plist)
    {
      ReadParameters_();
    }

    std::string variable;
    std::string region;
    std::string functional;
    const Teuchos::ParameterList& plist_;
  };

  // constructor and destructor
  Unstructured_observations(Teuchos::ParameterList& obs_list,
                            Amanzi::ObservationData& observation_data,
			    Epetra_MpiComm* comm);

  ~Unstructured_observations() {
    if (vo_ != NULL) delete vo_;
  }
  
  void RegisterComponentNames(std::vector<std::string> comp_names, int num_liquid) {
    comp_names_ = comp_names;
    num_liquid_ = num_liquid;
  }

  int MakeObservations(State& S);

  bool DumpRequested(const int);
  bool DumpRequested(const double);
  bool DumpRequested(const int, const double);
     
  void RegisterWithTimeStepManager(const Teuchos::Ptr<TimeStepManager>& tsm);

  void FlushObservations();

 private:
  double CalculateWaterTable_(State& S, AmanziMesh::Entity_ID_List& ids);

 protected:
  VerboseObject* vo_;
  
 private:
  int rank_;
  Teuchos::ParameterList obs_list_;
  Amanzi::ObservationData& observation_data_;
  std::map<std::string, Observable> observations;

  std::vector<std::string> comp_names_;
  int num_liquid_;
};

}  // namespace Amanzi

#endif
