#ifndef _UNSTRUCTURED_OBSERVATIONS_HPP_
#define _UNSTRUCTURED_OBSERVATIONS_HPP_

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Array.hpp"

#include "ObservationData.H"
#include "state.hh"

namespace Amanzi {

class UnstructuredObservations {

public:

  struct Observable {
    Observable(std::string state_id_,
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


  UnstructuredObservations(Teuchos::ParameterList observations_plist,
                           ObservationData& observation_data);

  void MakeObservations(const State& state);

private:

  Amanzi::ObservationData& observation_data_;

  std::map<std::string, Observable> observations_;

};

}


#endif
