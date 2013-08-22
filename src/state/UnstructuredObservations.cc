 /* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
Amanzi

License:
Author: Markus Berndt
        Ethan Coon (ecoon@lanl.gov)

Collection of Observations on an unstructured mesh.

------------------------------------------------------------------------- */

#include "UnstructuredObservations.hh"

#include "dbc.hh"
#include "errors.hh"
#include "exceptions.hh"
#include "State.hh"

#include <map>

namespace Amanzi {

UnstructuredObservations::UnstructuredObservations(
          Teuchos::ParameterList plist,
          const Teuchos::RCP<ObservationData>& observation_data,
          Epetra_MpiComm* comm) :
    observation_data_(observation_data) {

  // interpret paramerter list
  // loop over the sublists and create an observation for each
  for (Teuchos::ParameterList::ConstIterator i = plist.begin(); i != plist.end(); i++) {
    if (plist.isSublist(plist.name(i))) {
      Teuchos::ParameterList sublist = plist.sublist(plist.name(i));
      Teuchos::RCP<Observable> obs = Teuchos::rcp(new Observable(sublist, comm));

      observations_.insert(std::make_pair(sublist.name(), obs));
    }
  }
}

bool UnstructuredObservations::DumpRequested(int cycle, double time) const {
  for (ObservableMap::const_iterator lcv = observations_.begin();
       lcv != observations_.end(); ++lcv) {
    if (lcv->second->DumpRequested(cycle, time)) return true;
  }
  return false;
}

void UnstructuredObservations::MakeObservations(const State& S) {
  // loop over all observables
  for (ObservableMap::iterator lcv = observations_.begin();
       lcv != observations_.end(); ++lcv) {

    if (lcv->second->DumpRequested(S.cycle(), S.time())) {
      // data structure to store the observation
      Amanzi::ObservationData::DataTriple data_triplet;

      // make the observation
      lcv->second->Update(S, data_triplet);

      // push back into observation_data
      std::vector<Amanzi::ObservationData::DataTriple> &od =
          (*observation_data_)[lcv->first];
      od.push_back(data_triplet);
    }
  }
}

void UnstructuredObservations::RegisterWithTimeStepManager(
    const Teuchos::Ptr<TimeStepManager>& tsm) {

  // loop over all observables
  for (ObservableMap::iterator lcv = observations_.begin();
       lcv != observations_.end(); ++lcv) {
    // register
    lcv->second->RegisterWithTimeStepManager(tsm);
  }
}

}  // namespace Amanzi
