/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Markus Berndt
      Ethan Coon (coonet@ornl.gov)
*/

//! <MISSING_ONELINE_DOCSTRING>

#include "UnstructuredObservations.hh"

#include "State.hh"
#include "dbc.hh"
#include "errors.hh"
#include "exceptions.hh"

#include <map>

namespace Amanzi {

UnstructuredObservations::UnstructuredObservations(
    const Teuchos::RCP<Teuchos::ParameterList>& plist,
    const Comm_ptr_type& comm)
{
  // interpret paramerter list
  // loop over the sublists and create an observation for each
  for (auto& i : *plist) {
    if (plist->isSublist(i.first)) {
      auto sublist = Teuchos::sublist(plist, i.first);
      auto obs = Teuchos::rcp(new Observable(sublist, comm));
      observations_.insert(std::make_pair(sublist->name(), obs));
    }
  }
}

bool
UnstructuredObservations::DumpRequested(int cycle, double time) const
{
  for (ObservableMap::const_iterator lcv = observations_.begin();
       lcv != observations_.end();
       ++lcv) {
    if (lcv->second->DumpRequested(cycle, time)) return true;
  }
  return false;
}

void
UnstructuredObservations::MakeObservations(const State& S)
{
  // loop over all observables
  for (ObservableMap::iterator lcv = observations_.begin();
       lcv != observations_.end();
       ++lcv) {
    if (lcv->second->DumpRequested(S.cycle(), S.time())) {
      // data structure to store the observation
      Amanzi::ObservationData::DataQuadruple data_triplet;

      // make the observation
      lcv->second->Update(S, data_triplet);
    }
  }
}

void
UnstructuredObservations::RegisterWithTimeStepManager(
  TimeStepManager& tsm)
{
  // loop over all observables
  for (auto& lcv : observations_) {
    lcv.second->RegisterWithTimeStepManager(tsm);
  }
}

// It's not clear to me that this is necessary -- it seems that ofstream's
// destructor SHOULD flush (as in fclose), but maybe it doesn't?  Better safe
// than sorry...
void
UnstructuredObservations::Flush()
{
  for (ObservableMap::const_iterator lcv = observations_.begin();
       lcv != observations_.end();
       ++lcv) {
    lcv->second->Flush();
  }
}

} // namespace Amanzi
