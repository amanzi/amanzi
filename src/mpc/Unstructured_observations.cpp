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
  for (Teuchos::ParameterList::ConstIterator i = observations_plist_.begin(); i != observations_plist_.end(); i++) {

    if (observations_plist_.isSublist(observations_plist_.name(i))) {
      Teuchos::ParameterList observable_plist = observations_plist_.sublist(observations_plist_.name(i));

      std::vector<double> times;
      std::vector<std::vector<double> > time_sps;

      // get the observation times
      if (observable_plist.isSublist("time start period stop")) {
        Teuchos::ParameterList& tsps_list = observable_plist.sublist("time start period stop");
        for (Teuchos::ParameterList::ConstIterator it = tsps_list.begin(); it != tsps_list.end(); ++it) {
          std::string name = it->first;
          if (tsps_list.isSublist(name)) {
            Teuchos::ParameterList& itlist = tsps_list.sublist(name);
            if (itlist.isParameter("start period stop")) {
              Teuchos::Array<double> sps = itlist.get<Teuchos::Array<double> >("start period stop");
              time_sps.push_back(sps.toVector());
            }
          }
        }
      }
      if (observable_plist.isParameter("times")) {
        Teuchos::Array<double> vtimes = observable_plist.get<Teuchos::Array<double> >("times");
        times = vtimes.toVector();
      }

      std::vector<int> cycles;
      std::vector<std::vector<int> > cycle_sps;

      // get the observation cycles
      if (observable_plist.isSublist("cycle start period stop")) {
        Teuchos::ParameterList& csps_list = observable_plist.sublist("cycle start period stop");
        for (Teuchos::ParameterList::ConstIterator it = csps_list.begin(); it != csps_list.end(); ++it) {
          std::string name = it->first;
          if (csps_list.isSublist(name)) {
            Teuchos::ParameterList& itlist = csps_list.sublist(name);
            if (itlist.isParameter("start period stop")) {
              Teuchos::Array<int> csps = itlist.get<Teuchos::Array<int> >("start period stop");
              cycle_sps.push_back(csps.toVector());
            }
          }
        }
      }
      if (observable_plist.isParameter("cycles")) {
        Teuchos::Array<int> vcycles = observable_plist.get<Teuchos::Array<int> >("cycles");
        cycles = vcycles.toVector();
      }

      // loop over all variables listed and create an observable for each
      std::string var = observable_plist.get<string>("Variable");
      observations.insert(std::pair
                          <std::string, Observable>(observations_plist_.name(i),
                                                    Observable(var,
                                                               observable_plist.get<string>("Region"),
                                                               observable_plist.get<string>("Functional"),
                                                               cycles, cycle_sps, times, time_sps)));
    }
  }
}


void Unstructured_observations::make_observations(State& state)
{
  // loop over all observables
  for (std::map<std::string, Observable>::iterator i = observations.begin(); i != observations.end(); i++) {
    // for now we can only observe Integrals and Values

    if ( (i->second).functional != "Observation Data: Integral"  &&
         (i->second).functional != "Observation Data: Point" )  {
      Errors::Message m("Unstructured_observations: can only handle Functional == Observation Data: Integral, or Functional == Observation Data: Point");
      Exceptions::amanzi_throw(m);
    }

    std::string label = i->first;

    // make sure that we need to make an observation now
    if (observation_requested(state.get_time(), state.get_last_time(), (i->second).times, (i->second).sps) ||
        observation_requested(state.get_cycle(), (i->second).cycles, (i->second).csps ) ) {
      // we need to make an observation for each variable in the observable
      std::string var = (i->second).variable;

      // data structure to store the observation
      Amanzi::ObservationData::DataTriple data_triplet;

      // build the name of the observation
      std::stringstream ss;
      ss << label << ", " << var;

      std::vector<Amanzi::ObservationData::DataTriple> &od = observation_data[label];  //ss.str() ];

      if ((i->second).functional == "Observation Data: Integral")  {
        if (var == "Water") {
          data_triplet.value = state.water_mass();
        }
      } else if ((i->second).functional == "Observation Data: Point") {
        data_triplet.value = state.point_value((i->second).region, var);
      }

      data_triplet.is_valid = true;
      data_triplet.time = state.get_time();

      od.push_back(data_triplet);
    }
  }
}

void Unstructured_observations::register_with_time_step_manager(TimeStepManager& TSM) {
  // loop over all observables
  for (std::map<std::string, Observable>::const_iterator i = observations.begin();
       i != observations.end(); i++) {
    if ((i->second).sps.size() > 0) {
      for (std::vector<std::vector<double> >::const_iterator j=(i->second).sps.begin();
           j!=(i->second).sps.end(); ++j) {
        if (j->size() == 3) {
          TSM.RegisterTimeEvent((*j)[0], (*j)[1], (*j)[2]);
        }
      }
    }
    if ((i->second).times.size() > 0) {
      TSM.RegisterTimeEvent((i->second).times);
    }
  }
}

bool Unstructured_observations::observation_requested(double time, double last_time,
                                                      const std::vector<double>& T,
                                                      const std::vector<std::vector<double> >& SPS) {
  for (int i = 0; i < T.size(); i++)
    if (Amanzi::near_equal(T[i],time)) {
      return true;
    }
  if (SPS.size() > 0) {
    for (std::vector<std::vector<double> >::const_iterator i=SPS.begin(); i!=SPS.end(); ++i) {
      if  ( (Amanzi::near_equal(time,(*i)[0]) || time >= (*i)[0]) && 
	    (Amanzi::near_equal((*i)[2],-1.0) || time <= (*i)[2] || Amanzi::near_equal(time,(*i)[2]) ) ) {
        if (Amanzi::near_equal(time,(*i)[0])) {
	  return true;
	}
	double n_per_tmp = (time - (*i)[0])/(*i)[1];
        double n_periods = floor(n_per_tmp);
	if (Amanzi::near_equal(n_periods+1.0,n_per_tmp)) n_periods += 1.0;
        double tmp = (*i)[0] + n_periods*(*i)[1];
	if (Amanzi::near_equal(time,tmp)) {
	  return true;
	}
      }
    }
  }
  return false;
}


bool Unstructured_observations::observation_requested(int cycle,
                                                      const std::vector<int>& cyc,
                                                      const std::vector<std::vector<int> >& isps) {
  for (int i = 0; i < cyc.size(); i++)
    if (cyc[i] == cycle) return true;
  if (isps.size() > 0) {
    for (std::vector<std::vector<int> >::const_iterator i=isps.begin(); i!=isps.end(); ++i) {
      if  (cycle >= (*i)[0] && ((*i)[2] == -1 || cycle <= (*i)[2])) {
        if (  ( cycle-(*i)[0] ) % (*i)[1] == 0 ) return true;
      }
    }
  }

  return false;
}



}  // namespace Amanzi
