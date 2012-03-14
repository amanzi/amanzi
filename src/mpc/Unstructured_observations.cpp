#include "Unstructured_observations.hpp"

#include "dbc.hh"
#include "errors.hh"
#include "exceptions.hh"

#include <map>

namespace Amanzi {

  Unstructured_observations::Unstructured_observations(Teuchos::ParameterList observations_plist_,
						       Amanzi::ObservationData& observation_data_):
    observation_data(observation_data_)
  {
    // interpret paramerter list
    // loop over the sublists and create an observation for each
    for (Teuchos::ParameterList::ConstIterator i = observations_plist_.begin(); 
	 i != observations_plist_.end(); 
	 i++)
      {
	//  sublists
	if (observations_plist_.isSublist(observations_plist_.name(i))) 
	  {
	    Teuchos::ParameterList observable_plist = observations_plist_.sublist(observations_plist_.name(i));
	    
	    // get the observation times
	    Teuchos::Array<double> empty(0);
	    Teuchos::Array<double> sps = observable_plist.get<Teuchos::Array<double> >("Start_Period_Stop", empty);

	    // get the observation time values
	    Teuchos::Array<double> times = observable_plist.get<Teuchos::Array<double> >("Values", empty);
	    
	    // loop over all variables listed and create an observable for each
	    std::string var = observable_plist.get<string>("Variable");
	    observations.insert(std::pair
				<std::string, Observable>(observations_plist_.name(i), 
							  Observable(var,
								     observable_plist.get<string>("Region"),
								     observable_plist.get<string>("Functional"),
								     times, sps)));
	  }
      }    
    
  }
  
  
  
  void Unstructured_observations::make_observations(State& state)
  {

    // loop over all observables
    for (std::map<std::string, Observable>::iterator i = observations.begin();
	 i != observations.end();
	 i++) {
      
      // for now we can only observe Integrals and Values
      
      if ( (i->second).functional != "Observation Data: Integral"  &&  
	   (i->second).functional != "Observation Data: Point" )  {
	Errors::Message m("Unstructured_observations: can only handle Functional == Observation Data: Integral, or Functional == Observation Data: Point");
	Exceptions::amanzi_throw(m); 
      }
      
      
      std::string label = i->first;
      
      // make sure that we need to make an observation now
      if (observation_requested(state.get_time(), state.get_last_time(), (i->second).times, (i->second).sps)) {
	// we need to make an observation for each variable in the observable
	std::string var = (i->second).variable;
	
	// data structure to store the observation
	Amanzi::ObservationData::DataTriple data_triplet;
	
	// build the name of the observation
	std::stringstream ss;
	ss << label << ", " << var;
	
	std::vector<Amanzi::ObservationData::DataTriple> &od = observation_data[ label ]; //ss.str() ];
	
	if ((i->second).functional == "Observation Data: Integral")  {
	  if (var == "Water") {
	    data_triplet.value   = state.water_mass();
	  }
	} else if ((i->second).functional == "Observation Data: Point") {
	  data_triplet.value   = state.point_value((i->second).region, var);
	}
	
	data_triplet.is_valid = true;
	data_triplet.time = state.get_time();
	
	od.push_back(data_triplet);
      }
      
    }
  } 
  
  bool  Unstructured_observations::observation_requested(double time, double last_time, Teuchos::Array<double>& T, Teuchos::Array<double>& SPS) {
    
    for (int i=0; i<T.size(); i++) {
      if ( last_time < T[i] && T[i] <= time ) {
	return true;
      }
    }
    
    if (SPS.size() > 0) {
      if ( time >= SPS[0] && (SPS[2]==-1.0 || time <=SPS[2]) ) {
	if (time == SPS[0]) return true;
	
	int n0 =  floor( (last_time - SPS[0])/SPS[1] );
	int n1 =  floor( (time - SPS[0])/SPS[1] );
	
	if (n0+1 == n1) return true;
      }
    }
    
    return false;  
  }
  


}
