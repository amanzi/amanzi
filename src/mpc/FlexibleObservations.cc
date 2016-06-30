/*
  This is the flow component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Markus Berndt  
*/

#include <map>

// Amanzi
#include "dbc.hh"
#include "errors.hh"
#include "exceptions.hh"
#include "Point.hh"
#include "RegionPlane.hh"
#include "RegionPolygon.hh"
#include "ReconstructionCell.hh"
#include "Units.hh"

// MPC
#include "FlexibleObservations.hh"
#include "ObservableFactory.hh"
#include "Observable.hh"

namespace Amanzi {

/* ******************************************************************
* Constructor.
****************************************************************** */
FlexibleObservations::FlexibleObservations(
    Teuchos::RCP<Teuchos::ParameterList> coordinator_list,
    Teuchos::RCP<Teuchos::ParameterList> obs_list,
    Teuchos::RCP<Teuchos::ParameterList> units_list,
    Amanzi::ObservationData& observation_data,
    Teuchos::RCP<AmanziMesh::Mesh> mesh)
    : observation_data_(observation_data),
      obs_list_(obs_list),
      coordinator_list_(coordinator_list)
{
  rank_ = mesh->get_comm()->MyPID();

  // initialize units
  units_.Init(*units_list);

  Teuchos::ParameterList tmp_list;
  tmp_list.set<std::string>("verbosity level", "high");
  vo_ = new VerboseObject("Observations", tmp_list);

  // loop over the sublists and create an observation for each
  for (Teuchos::ParameterList::ConstIterator i = obs_list_->begin(); i != obs_list_->end(); i++) {

    if (obs_list_->isSublist(obs_list_->name(i))) {
      Teuchos::ParameterList observable_plist = obs_list_->sublist(obs_list_->name(i));

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
      std::string var = observable_plist.get<std::string>("variable");
      observations.insert(std::pair<std::string, Teuchos::RCP<Observable> >(
          obs_list_->name(i), 
	  CreateObservable(*coordinator_list,  observable_plist, *units_list, mesh)));
          // Observable(var, observable_plist.get<std::string>("region"),
          //            observable_plist.get<std::string>("functional"),
          //            observable_plist, comm)));
    }
  }
}


/* ******************************************************************
* Process data to extract observations.
****************************************************************** */
int FlexibleObservations::MakeObservations(State& S)
{
  Errors::Message msg;
  int num_obs(0);
  int dim = S.GetMesh()->space_dimension();

  // loop over all observables
  for (std::map<std::string, Teuchos::RCP<Observable> >::iterator i = observations.begin(); i != observations.end(); i++) {
    if ((i->second)->DumpRequested(S.time()) || (i->second)->DumpRequested(S.cycle())) {

      num_obs++;          
      
      std::string label = i->first;
      
      //we need to make an observation for each variable in the observable
      std::string var = (i->second)->variable_;
      
      // data structure to store the observation
      Amanzi::ObservationData::DataTriple data_triplet;
      
      // build the name of the observation
      std::stringstream ss;
      ss << label << ", " << var;
      
      std::vector<Amanzi::ObservationData::DataTriple>& od = observation_data_[label]; 

      double value(0.0), volume(0.0);

      (i->second) -> ComputeObservation(S, &value, &volume);

      if (var == "drawdown" || var=="permeability-weighted drawdown"){
        if (od.size() > 0) { 
          value = od.begin()->value * volume - value;
        }
      }
    
      // syncronize the result across processors
      double result;
      S.GetMesh()->get_comm()->SumAll(&value, &result, 1);
      
      double vresult;
      S.GetMesh()->get_comm()->SumAll(&volume, &vresult, 1);
 
      if ((i->second)->functional_ == "observation data: integral") {  
        data_triplet.value = result;  
      } else if ((i->second)->functional_ == "observation data: point") {
        data_triplet.value = result / vresult;
      }
      
      data_triplet.is_valid = true;
      data_triplet.time = S.time();

      bool time_exist = false;
      for (std::vector<Amanzi::ObservationData::DataTriple>::iterator it = od.begin(); it != od.end(); ++it) {
        if (it->time == data_triplet.time) {
          time_exist = true;
          break;
        }
      }
            
      if (!time_exist) od.push_back(data_triplet);
    }
  }

  FlushObservations();
  return num_obs;
}


/* ******************************************************************
* Save observation based on time or cycle.
****************************************************************** */
bool FlexibleObservations::DumpRequested(const double time) {
  bool result = false;
  for (std::map<std::string, Teuchos::RCP<Observable> >::iterator i = observations.begin(); i != observations.end(); i++) {
    result = result || (i->second)->DumpRequested(time);
  }  
  return result;
}


bool FlexibleObservations::DumpRequested(const int cycle) {
  bool result = false;
  for (std::map<std::string, Teuchos::RCP<Observable> >::iterator i = observations.begin(); i != observations.end(); i++) {
    result = result || (i->second)->DumpRequested(cycle);
  }  
  return result;  
}


bool FlexibleObservations::DumpRequested(const int cycle, const double time) {
  return DumpRequested(time) || DumpRequested(cycle);
}


/******************************************************************
* Loop over all observations and register each of them with the time 
* step manager.
******************************************************************/
void FlexibleObservations::RegisterWithTimeStepManager(const Teuchos::Ptr<TimeStepManager>& tsm) {
  for (std::map<std::string, Teuchos::RCP<Observable> >::iterator i = observations.begin(); i != observations.end(); i++) {
    (i->second)->RegisterWithTimeStepManager(tsm);
  }  
}


/* ******************************************************************
* Write observatins to a file. Clsoe the file to flush data.
****************************************************************** */
void FlexibleObservations::FlushObservations()
{
  if (obs_list_->isParameter("observation output filename")) {
    std::string obs_file = obs_list_->get<std::string>("observation output filename");
    int precision = obs_list_->get<int>("precision", 16);

    if (rank_ == 0) {
      std::ofstream out;
      out.open(obs_file.c_str(), std::ios::out);
      
      out.precision(precision);
      out.setf(std::ios::scientific);

      out << "Observation Name, Region, Functional, Variable, Time, Value\n";
      out << "===========================================================\n";

      for (Teuchos::ParameterList::ConstIterator i = obs_list_->begin(); i != obs_list_->end(); ++i) {
        std::string label = obs_list_->name(i);
        const Teuchos::ParameterEntry& entry = obs_list_->getEntry(label);
        if (entry.isList()) {
          const Teuchos::ParameterList& ind_obs_list = obs_list_->sublist(label);
          std::vector<Amanzi::ObservationData::DataTriple>& od = observation_data_[label]; 

          for (int j = 0; j < od.size(); j++) {
            if (od[j].is_valid) {
              if (!out.good()) {
                std::cout << "PROBLEM BEFORE" << std::endl;
              }
              std::string var = ind_obs_list.get<std::string>("variable");
              out << label << ", "
                  << ind_obs_list.get<std::string>("region") << ", "
                  << ind_obs_list.get<std::string>("functional") << ", "
                  << var << ", "
                  << od[j].time << ", "
                  << (((var == "permeability-weighted drawdown" || 
                        var == "drawdown") && !j) ? 0.0 : od[j].value) << '\n';
              if (!out.good()) {
                std::cout << "PROBLEM AFTER" << std::endl;
              }
            }
          }
        }
      }
      out.close();
    }
  }
}

}  // namespace Amanzi
