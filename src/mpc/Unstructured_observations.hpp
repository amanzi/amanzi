#ifndef _UNSTRUCTURED_OBSERVATIONS_HPP_
#define _UNSTRUCTURED_OBSERVATIONS_HPP_

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Array.hpp"

#include "ObservationData.H"
#include "State.hpp"
#include "time_step_manager.hh"

namespace Amanzi {


class Unstructured_observations
{

 public:

  struct Observable
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
        sps(sps_), cycles(cycles_), csps(csps_)
    {};

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

  void make_observations(State& state);
  bool observation_requested(double time, double last_time, const std::vector<double>& T,
                             const std::vector<std::vector<double> >& SPS);
  bool observation_requested(int cycle, const std::vector<int>& cyc, 
			     const std::vector<std::vector<int> >& csps);
  void register_with_time_step_manager(TimeStepManager& TSM);
 private:

  Amanzi::ObservationData& observation_data;

  std::map<std::string, Observable> observations;

};

}


#endif
