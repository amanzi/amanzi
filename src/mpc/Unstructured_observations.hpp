#ifndef _UNSTRUCTURED_OBSERVATIONS_HPP_
#define _UNSTRUCTURED_OBSERVATIONS_HPP_

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Array.hpp"

#include "ObservationData.H"
#include "State.hpp"

namespace Amanzi {


  class Unstructured_observations
  {
    
  public:

    struct Observable 
    {
      Observable (std::string variable_,
		  std::string region_,
		  std::string functional_,
		  Teuchos::Array<double> times_,
		  Teuchos::Array<double> sps_):
	variable(variable_), region(region_),
	functional(functional_), times(times_),
	sps(sps_)
      {};
      
      std::string variable;
      std::string region;
      std::string functional;
      Teuchos::Array<double> times;
      Teuchos::Array<double> sps;   // start period stop
    };


    Unstructured_observations (Teuchos::ParameterList observations_plist_,
			       Amanzi::ObservationData& observation_data_);

    void make_observations(State& state);
    bool observation_requested(double time, double last_time, Teuchos::Array<double>& T, Teuchos::Array<double>& SPS);
  private:
    
    Amanzi::ObservationData& observation_data;

    std::map<std::string, Observable> observations;

  };

}
			       

#endif
