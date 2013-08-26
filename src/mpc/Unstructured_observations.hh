#ifndef _UNSTRUCTURED_OBSERVATIONS_HPP_
#define _UNSTRUCTURED_OBSERVATIONS_HPP_

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Array.hpp"

#include "ObservationData.hh"
#include "State.hh"
#include "io_event.hh"

namespace Amanzi {


class Unstructured_observations {

 public:

  struct Observable : public IOEvent
  {
    Observable (std::string variable_,
                std::string region_,
                std::string functional_,
                std::vector<int> cycles_,
                std::vector<std::vector<int> > csps_,
                std::vector<double> times_,
                std::vector<std::vector<double> > sps_):
        variable(variable_), region(region_),
        functional(functional_), times(times_),
        sps(sps_), cycles(cycles_), csps(csps_) {}

    std::string variable;
    std::string region;
    std::string functional;
    std::vector<double> times;
    std::vector<int> cycles;
    std::vector<std::vector<double> > sps;   // time start period stop
    std::vector<std::vector<int> > csps;  // cycle start period stop
  };


  Unstructured_observations (Teuchos::ParameterList observations_plist_,
                             Amanzi::ObservationData& observation_data_);
  
  void register_component_names(std::vector<std::string> comp_names) {
    comp_names_ = comp_names;
  }

  void make_observations(State& state);

  bool DumpRequested(const int);
  bool DumpRequested(const double);
  bool DumpRequested(const int, const double);
     
  void RegisterWithTimeStepManager(const Teuchos::Ptr<TimeStepManager>& tsm);
  

 private:

  Amanzi::ObservationData& observation_data;
  std::map<std::string, Observable> observations;
  std::vector<std::string> comp_names_;

};

}


#endif
