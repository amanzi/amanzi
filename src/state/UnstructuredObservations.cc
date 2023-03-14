/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Markus Berndt
           Ethan Coon (ecoon@lanl.gov)
*/

/*
  This is the state component of the Amanzi code.

  Collection of Observations on an unstructured mesh.
*/

#include <map>
#include <memory>
#include "boost/filesystem/operations.hpp"

#include "UnstructuredObservations.hh"

#include "dbc.hh"
#include "errors.hh"
#include "exceptions.hh"
#include "State.hh"

namespace Amanzi {

UnstructuredObservations::UnstructuredObservations(Teuchos::ParameterList& plist)
  : IOEvent(plist),
    write_(false),
    num_total_(0),
    count_(0),
    time_integrated_(false),
    observed_once_(false)
{
  // file information
  filename_ = plist.get<std::string>("observation output filename");
  delimiter_ = plist.get<std::string>("delimiter", ",");
  interval_ = plist.get<int>("write interval", 1);
  time_unit_ = plist.get<std::string>("time units", "s");
  Utils::Units unit;
  bool flag = false;
  time_unit_factor_ = unit.ConvertTime(1., "s", time_unit_, flag);

  // this sets the communicator
  writing_domain_ = plist.get<std::string>("domain", "domain");

  // interpret parameter list
  // loop over the sublists and create an observable for each
  if (plist.isSublist("observed quantities")) {
    Teuchos::ParameterList& oq_list = plist.sublist("observed quantities");
    for (auto it : oq_list) {
      if (oq_list.isSublist(it.first)) {
        auto obs = Teuchos::rcp(new Observable(oq_list.sublist(it.first)));
        observables_.emplace_back(obs);
        time_integrated_ |= obs->is_time_integrated();
      } else {
        Errors::Message msg;
        msg << "Observation list \"observed quantities\" should contain only sublists.";
        Exceptions::amanzi_throw(msg);
      }
    }

    if (observables_.size() == 0) {
      Errors::Message msg;
      msg << "Observation list \"observed quantities\" is empty.";
      Exceptions::amanzi_throw(msg);
    }

  } else {
    // old style, single list/single entry
    auto obs = Teuchos::rcp(new Observable(plist));
    observables_.emplace_back(obs);
    time_integrated_ |= obs->is_time_integrated();
  }
}

void
UnstructuredObservations::Setup(const Teuchos::Ptr<State>& S)
{
  // set the communicator
  comm_ = Teuchos::null;
  if (S->HasMesh(writing_domain_)) comm_ = S->GetMesh(writing_domain_)->get_comm();

  // require fields, evaluators
  for (auto& obs : observables_) {
    obs->set_comm(comm_);
    obs->Setup(S);
  }

  // what rank writes the file?
  write_ = false;
  if (comm_ != Teuchos::null && comm_->MyPID() == 0) write_ = true;
}

void
UnstructuredObservations::MakeObservations(const Teuchos::Ptr<State>& S)
{
  if (comm_ == Teuchos::null) return;

  if (!observed_once_) {
    // final setup, open file handle, etc
    for (auto& obs : observables_) obs->FinalizeStructure(S);

    num_total_ = 0;
    for (const auto& obs : observables_) {
      int num_local = obs->get_num_vectors();
      int num_global = -1;
      comm_->MaxAll(&num_local, &num_global, 1);
      num_total_ += num_global;
      if (num_global <= 0) {
        // at last we can finally catch a bad observation name here -- no one
        // has this variable, it does not exist!
        Errors::Message msg;
        msg << "Invalid observation requested for variable \"" << obs->get_variable()
            << "\"; no such field exists.";
        Exceptions::amanzi_throw(msg);
      }
    }
    integrated_observation_.resize(num_total_);
    observed_once_ = true;

    if (write_) InitFile_(S);
  }

  bool dump_requested = DumpRequested(S->get_cycle(), S->get_time());
  if (time_integrated_) {
    if (dump_requested) {
      std::vector<double> observation(num_total_, Observable::nan);
      int loc_start = 0;
      for (auto& obs : observables_) {
        obs->Update(S, observation, loc_start);
        if (obs->is_time_integrated()) {
          for (int i = 0; i != obs->get_num_vectors(); ++i) {
            observation[loc_start + i] += integrated_observation_[loc_start + i];
            integrated_observation_[loc_start + i] = 0.;
          }
        }
        loc_start += obs->get_num_vectors();
      }

      // write
      if (write_) Write_(S->get_time() * time_unit_factor_, observation);
    } else {
      std::vector<double> observation(num_total_, Observable::nan);

      // loop over all observables
      int loc_start = 0;
      for (auto& obs : observables_) {
        if (obs->is_time_integrated()) {
          obs->Update(S, observation, loc_start);
          for (int i = 0; i != obs->get_num_vectors(); ++i) {
            integrated_observation_[loc_start + i] += observation[loc_start + i];
          }
        }
        loc_start += obs->get_num_vectors();
      }
    }
  } else if (dump_requested) {
    std::vector<double> observation(num_total_, Observable::nan);

    // loop over all observables
    int loc_start = 0;
    for (auto& obs : observables_) {
      obs->Update(S, observation, loc_start);
      loc_start += obs->get_num_vectors();
    }

    // write
    if (write_) Write_(S->get_time() * time_unit_factor_, observation);
  }
}


//
// Note this should only be called on one rank
//
// Note also that this (along with Write_) may become a separate class for
// different file types (e.g. text vs netcdf)
//
void
UnstructuredObservations::InitFile_(const Teuchos::Ptr<const State>& S)
{
  if (boost::filesystem::portable_file_name(filename_)) {
    fid_ = std::make_unique<std::ofstream>(filename_.c_str());

    *fid_ << "# Observation File: " << filename_ << " column names:" << std::endl
          << "# -----------------------------------------------------------------------------"
          << std::endl
          << "# Observation Name: time [" << time_unit_ << "]" << std::endl;

    for (const auto& obs : observables_) {
      *fid_ << "# -----------------------------------------------------------------------------"
            << std::endl
            << "# Observation Name: " << obs->get_name() << std::endl
            << "# Region: " << obs->get_region() << std::endl
            << "# Reduction: " << obs->get_reduction() << std::endl;
      if (!obs->get_modifier().empty()) {
        *fid_ << "# Modifier:" << std::endl
              << "#   " << obs->get_modifier() << std::endl;
      }
      *fid_
            << "# Variable: " << obs->get_variable() << std::endl
            << "# Number of Vectors: " << obs->get_num_vectors() << std::endl;

      if (obs->get_degree_of_freedom() >= 0)
        *fid_ << "# DoF: " << obs->get_degree_of_freedom() << std::endl;
    }
    *fid_ << "# ============================================================================="
          << std::endl;
    *fid_ << "\"time [" << time_unit_ << "]\"";
    for (const auto& obs : observables_) {
      if (obs->get_num_vectors() > 1) {
        auto subfield_names = S->GetRecordSet(obs->get_variable()).subfieldnames();
        std::vector<std::string> default_subfieldnames;
        if (!subfield_names) {
          for (int i = 0; i != obs->get_num_vectors(); ++i)
            default_subfieldnames.emplace_back(std::string("dof ") + std::to_string(i));
          subfield_names = &default_subfieldnames;
        }

        if (obs->get_degree_of_freedom() >= 0) {
          *fid_ << delimiter_ << "\"" << obs->get_name() << " "
                << (*subfield_names)[obs->get_degree_of_freedom()] << "\"";
        } else {
          for (int i = 0; i != obs->get_num_vectors(); ++i) {
            *fid_ << delimiter_ << "\"" << obs->get_name() << " " << (*subfield_names)[i] << "\"";
          }
        }
      } else {
        *fid_ << delimiter_ << "\"" << obs->get_name() << "\"";
      }
    }
    *fid_ << std::endl << std::scientific;
    fid_->precision(12);
  } else {
    Errors::Message msg;
    msg << "Invalid filename for observation: \"" << filename_ << "\"";
    Exceptions::amanzi_throw(msg);
  }
}


void
UnstructuredObservations::Write_(double time, const std::vector<double>& obs)
{
  if (fid_.get()) {
    *fid_ << time;
    for (auto val : obs) *fid_ << delimiter_ << val;
    *fid_ << std::endl;
    ++count_;
    if (count_ % interval_ == 0) fid_->flush();
  }
}


// It's not clear to me that this is necessary -- it seems that ofstream's
// destructor SHOULD flush (as in fclose), but maybe it doesn't?  Better safe
// than sorry...
void
UnstructuredObservations::Flush()
{
  if (fid_.get()) fid_->flush();
}

} // namespace Amanzi
