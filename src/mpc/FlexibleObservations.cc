/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Markus Berndt (berndt@lanl.gov)
           Daniil Svyatskiy
*/

/*
  Multi-Process Coordinator

*/

#include <map>

// Amanzi
#include "dbc.hh"
#include "errors.hh"
#include "exceptions.hh"
#include "Point.hh"
#include "RegionPlane.hh"
#include "RegionPolygon.hh"
#include "Units.hh"

// MPC
#include "FlexibleObservations.hh"
#include "ObservableFactory.hh"
#include "ObservableAmanzi.hh"

namespace Amanzi {

/* ******************************************************************
* Constructor.
****************************************************************** */
FlexibleObservations::FlexibleObservations(Teuchos::RCP<Teuchos::ParameterList> coordinator_list,
                                           Teuchos::RCP<Teuchos::ParameterList> obs_list,
                                           Teuchos::RCP<Teuchos::ParameterList> units_list,
                                           Amanzi::ObservationData& observation_data,
                                           Teuchos::RCP<const State> S)
  : obs_list_(obs_list), coordinator_list_(coordinator_list), observation_data_(observation_data)
{
  rank_ = S->GetMesh("domain")->get_comm()->MyPID();

  // initialize units
  units_.Init(*units_list);

  Teuchos::ParameterList tmp_list;
  tmp_list.set<std::string>("verbosity level", "high");
  vo_ = new VerboseObject("Observations", tmp_list);

  // loop over the sublists and create an observation for each
  for (auto i = obs_list_->begin(); i != obs_list_->end(); i++) {
    if (obs_list_->isSublist(obs_list_->name(i))) {
      Teuchos::ParameterList observable_plist = obs_list_->sublist(obs_list_->name(i));

      std::vector<double> times;
      std::vector<std::vector<double>> time_sps;

      // get the observation times
      if (observable_plist.isSublist("time start period stop")) {
        Teuchos::ParameterList& tsps_list = observable_plist.sublist("time start period stop");
        for (auto it = tsps_list.begin(); it != tsps_list.end(); ++it) {
          std::string name = it->first;
          if (tsps_list.isSublist(name)) {
            Teuchos::ParameterList& itlist = tsps_list.sublist(name);
            if (itlist.isParameter("start period stop")) {
              Teuchos::Array<double> sps = itlist.get<Teuchos::Array<double>>("start period stop");
              time_sps.push_back(sps.toVector());
            }
          }
        }
      }
      if (observable_plist.isParameter("times")) {
        Teuchos::Array<double> vtimes = observable_plist.get<Teuchos::Array<double>>("times");
        times = vtimes.toVector();
      }

      std::vector<int> cycles;
      std::vector<std::vector<int>> cycle_sps;

      // get the observation cycles
      if (observable_plist.isSublist("cycle start period stop")) {
        Teuchos::ParameterList& csps_list = observable_plist.sublist("cycle start period stop");
        for (auto it = csps_list.begin(); it != csps_list.end(); ++it) {
          std::string name = it->first;
          if (csps_list.isSublist(name)) {
            Teuchos::ParameterList& itlist = csps_list.sublist(name);
            if (itlist.isParameter("start period stop")) {
              Teuchos::Array<int> csps = itlist.get<Teuchos::Array<int>>("start period stop");
              cycle_sps.push_back(csps.toVector());
            }
          }
        }
      }

      if (observable_plist.isParameter("cycles")) {
        Teuchos::Array<int> vcycles = observable_plist.get<Teuchos::Array<int>>("cycles");
        cycles = vcycles.toVector();
      }

      // loop over all variables listed and create an observable for each
      std::string var = observable_plist.get<std::string>("variable");
      Key domain_name = observable_plist.get<std::string>("domain name", "domain");
      observations.insert(std::pair<std::string, Teuchos::RCP<Observable>>(
        obs_list_->name(i),
        CreateObservable(
          *coordinator_list, observable_plist, *units_list, S->GetMesh(domain_name))));
      // Observable(var, observable_plist.get<std::string>("region"),
      //            observable_plist.get<std::string>("functional"),
      //            observable_plist, comm)));
    }
  }
}


/* ******************************************************************
* Process data to extract observations.
****************************************************************** */
int
FlexibleObservations::MakeObservations(State& S)
{
  int num_obs(0);
  std::string unit;

  // loop over all observables
  for (std::map<std::string, Teuchos::RCP<Observable>>::iterator i = observations.begin();
       i != observations.end();
       i++) {
    if ((i->second)->DumpRequested(S.get_time()) || (i->second)->DumpRequested(S.get_cycle())) {
      num_obs++;

      // we need to make an observation for each variable in the observable
      std::string var = (i->second)->variable_;

      // data structure to store the observation
      Amanzi::ObservationData::DataQuadruple data_quad;

      std::string label = i->first;
      std::vector<Amanzi::ObservationData::DataQuadruple>& od = observation_data_[label];

      double value(0.0), volume(0.0);
      i->second->ComputeObservation(S, &value, &volume, unit, 0.0);

      if (var == "drawdown" || var == "permeability-weighted drawdown") {
        if (od.size() > 0) { value = od.begin()->value * volume - value; }
      }

      // syncronize the result across processors
      double data_out[2], data_in[2] = { value, volume };
      S.GetMesh()->get_comm()->SumAll(data_in, data_out, 2);

      if ((i->second)->functional_ == "observation data: integral") {
        data_quad.value = data_out[0];
        unit.append("*m^3");
      } else if ((i->second)->functional_ == "observation data: point") {
        data_quad.value = data_out[0] / std::max(1e-16, data_out[1]);
      }

      data_quad.is_valid = true;
      data_quad.time = S.get_time();
      data_quad.unit = unit;

      bool time_exist = false;
      for (std::vector<Amanzi::ObservationData::DataQuadruple>::iterator it = od.begin();
           it != od.end();
           ++it) {
        if (it->time == data_quad.time) {
          time_exist = true;
          break;
        }
      }

      if (!time_exist) od.push_back(data_quad);
    }
  }

  FlushObservations();
  return num_obs;
}


/* ******************************************************************
* Process data to extract observations.
****************************************************************** */
int
FlexibleObservations::MakeContinuousObservations(State& S)
{
  for (auto i = observations.begin(); i != observations.end(); i++) {
    std::string var = i->second->variable_;
    if (var.find(" breakthrough curve") != std::string::npos) {
      double value(0.0), volume(0.0);
      std::string unit;
      double dt = S.final_time() - S.initial_time();
      i->second->ComputeObservation(S, &value, &volume, unit, dt);
    }
  }
  return 0;
}


/* ******************************************************************
* Save observation based on time or cycle.
****************************************************************** */
bool
FlexibleObservations::DumpRequested(const double time)
{
  bool result = false;
  for (std::map<std::string, Teuchos::RCP<Observable>>::iterator i = observations.begin();
       i != observations.end();
       i++) {
    result = result || (i->second)->DumpRequested(time);
  }
  return result;
}


bool
FlexibleObservations::DumpRequested(const int cycle)
{
  bool result = false;
  for (std::map<std::string, Teuchos::RCP<Observable>>::iterator i = observations.begin();
       i != observations.end();
       i++) {
    result = result || (i->second)->DumpRequested(cycle);
  }
  return result;
}


bool
FlexibleObservations::DumpRequested(const int cycle, const double time)
{
  return DumpRequested(time) || DumpRequested(cycle);
}


/******************************************************************
* Loop over all observations and register each of them with the time
* step manager.
******************************************************************/
void
FlexibleObservations::RegisterWithTimeStepManager(const Teuchos::Ptr<TimeStepManager>& tsm)
{
  for (std::map<std::string, Teuchos::RCP<Observable>>::iterator i = observations.begin();
       i != observations.end();
       i++) {
    (i->second)->RegisterWithTimeStepManager(tsm);
  }
}


/* ******************************************************************
* Write observatins to a file. Clsoe the file to flush data.
****************************************************************** */
void
FlexibleObservations::FlushObservations()
{
  bool flag1, flag2(true);

  if (obs_list_->isParameter("observation output filename")) {
    std::string obs_file = obs_list_->get<std::string>("observation output filename");
    int precision = obs_list_->get<int>("precision", 16);

    Utils::UnitsSystem system = units_.system();
    system.time = obs_list_->get<std::string>("time unit", system.time);
    system.mass = obs_list_->get<std::string>("mass unit", system.mass);
    system.length = obs_list_->get<std::string>("length unit", system.length);
    system.concentration = obs_list_->get<std::string>("concentration unit", system.concentration);

    if (rank_ == 0) {
      std::ofstream out;
      out.open(obs_file.c_str(), std::ios::out);

      out.precision(precision);
      out.setf(std::ios::scientific);

      out << "Observation Name, Region, Functional, Variable, Time [" << system.time << "], Value ["
          << system.mass << " " << system.length << " " << system.concentration << "]\n";
      out << "============================================================================\n";

      for (auto it = obs_list_->begin(); it != obs_list_->end(); ++it) {
        std::string label = obs_list_->name(it);
        const Teuchos::ParameterEntry& entry = obs_list_->getEntry(label);
        if (entry.isList()) {
          const Teuchos::ParameterList& ind_obs_list = obs_list_->sublist(label);
          std::vector<Amanzi::ObservationData::DataQuadruple>& od = observation_data_[label];

          for (int j = 0; j < od.size(); j++) {
            if (od[j].is_valid) {
              std::string var, name, out_unit;
              double out_value, mol_mass(1.0);

              // ugly way to extract molar mass of a specie
              var = ind_obs_list.get<std::string>("variable");
              std::stringstream ss(var);
              ss >> name;
              for (int i = 0; i < comp_names_.size(); ++i) {
                if (comp_names_[i] == name) {
                  mol_mass = comp_mol_masses_[i];
                  break;
                }
              }

              out_unit = units_.ConvertUnitS(od[j].unit, system);
              out_value = units_.ConvertUnitD(od[j].value, od[j].unit, out_unit, mol_mass, flag1);
              flag2 &= flag1;

              out << label << ", " << ind_obs_list.get<std::string>("region") << ", "
                  << ind_obs_list.get<std::string>("functional") << ", " << var << ", "
                  << units_.ConvertTime(od[j].time, "s", system.time, flag1) << ", "
                  << (((var == "permeability-weighted drawdown" || var == "drawdown") && !j) ?
                        0.0 :
                        out_value)
                  << '\n';
              flag2 &= flag1;
            }
          }
        }
      }
      out.close();
    }
  }

  if (!flag2) {
    Errors::Message msg;
    msg << "Conversion of units in observations has failed.\n";
    Exceptions::amanzi_throw(msg);
  }
}

} // namespace Amanzi
