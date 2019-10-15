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
  Teuchos::ParameterList& plist,
  const Teuchos::RCP<ObservationData>& observation_data, Epetra_MpiComm* comm)
  : observation_data_(observation_data)
{
  // interpret paramerter list
  // loop over the sublists and create an observation for each
  for (Teuchos::ParameterList::ConstIterator i = plist.begin();
       i != plist.end();
       i++) {
    if (plist.isSublist(plist.name(i))) {
      Teuchos::ParameterList sublist = plist.sublist(plist.name(i));
      Teuchos::RCP<Observable> obs =
        Teuchos::rcp(new Observable(sublist, comm));

      observations_.insert(std::make_pair(sublist.name(), obs));
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

      // push back into observation_data
      if (observation_data_ != Teuchos::null) {
        std::vector<Amanzi::ObservationData::DataQuadruple>& od =
          (*observation_data_)[lcv->first];
        od.push_back(data_triplet);
      }
    }
  }
}

void
UnstructuredObservations::RegisterWithTimeStepManager(
  const Teuchos::Ptr<TimeStepManager>& tsm)
{
  // loop over all observables
  for (ObservableMap::iterator lcv = observations_.begin();
       lcv != observations_.end();
       ++lcv) {
    // register
    lcv->second->RegisterWithTimeStepManager(tsm);
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
