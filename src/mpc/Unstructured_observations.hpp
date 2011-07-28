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
      Observable (std::string state_id_,
		  std::string region_,
		  std::string functional_,
		  Teuchos::Array<double> times_):
	state_id(state_id_), region(region_),
	functional(functional_), times(times_) {};

      std::string state_id;
      std::string region;
      std::string functional;
      Teuchos::Array<double> times;
    };


    Unstructured_observations (Teuchos::ParameterList observations_plist_,
			       Amanzi::ObservationData& observation_data_);

    void make_observations(State& state);

  private:
    
    Amanzi::ObservationData& observation_data;

    std::map<std::string, Observable> observations;

  };

}
			       

#endif
