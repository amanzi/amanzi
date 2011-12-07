#include "Unstructured_observations.hh"

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
	    Teuchos::Array<double> times = observable_plist.get<Teuchos::Array<double> >("times");

	    // create observation data structure
	    for (Teuchos::Array<double>::const_iterator itimes = times.begin(); 
		 itimes != times.end(); itimes++)
	      {
		Amanzi::ObservationData::DataTriple dt;
		dt.time = *itimes;
		observation_data[observations_plist_.name(i)].push_back(dt);
	      }
	    
	    observations.insert(std::pair
				<std::string,Observable>(observations_plist_.name(i),
							 Observable(observable_plist.get<string>("state id"),
								    observable_plist.get<string>("region"),
								    observable_plist.get<string>("functional"),
								    times)));
	  }
	else
	  {
	    Errors::Message m("Unstructured_observations: the Observation sublist contains an entry that is not a sublist!");
	    Exceptions::amanzi_throw(m);
	  }
      }    

  }
  
  
  
  void Unstructured_observations::make_observations(State& state)
  {
    for (std::map<std::string, Observable>::iterator i = observations.begin();
	 i != observations.end();
	 i++)

      {

	if ( (i->second).region != "all" )
	  {
	    Errors::Message m("Unstructured_observations: can only handle region == all");
	    Exceptions::amanzi_throw(m);
	  }

	if ( (i->second).state_id != "water" )
	  {
	    Errors::Message m("Unstructured_observations: can only handle state id == water");
	    Exceptions::amanzi_throw(m);
	  }
	    
	if ( (i->second).functional != "integral" )
	  {
	    Errors::Message m("Unstructured_observations: can only handle functional == integral");
	    Exceptions::amanzi_throw(m);
	  }
	
	
	// if ( (i->second).region == "all"          && 
	//      (i->second).state_id == "water"      &&
	//      (i->second).functional == "integral" ) 
	//   {
	//     std::string label = i->first;

	//     std::vector<Amanzi::ObservationData::DataTriple>::iterator it; 
	//     std::vector<Amanzi::ObservationData::DataTriple> &od = observation_data[label];
	    
	//     for ( it = od.begin(); it != od.end(); it++)
	//       {
	// 	if  ( state.get_time() >= it->time ) 
	//      	  {
	// 	    if ( ! it->is_valid ) 
	// 	      {
	//      		it->value   = state.water_mass();
	// 		it->is_valid = true;
	// 	      }
	//      	  }
	//       }
	//   }
       }
      }

}

