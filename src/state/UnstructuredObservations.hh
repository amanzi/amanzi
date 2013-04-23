 /* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
Amanzi

License:
Author: Markus Berndt
        Ethan Coon (ecoon@lanl.gov)

Collection of Observations on an unstructured mesh.

------------------------------------------------------------------------- */


#ifndef AMANZI_UNSTRUCTURED_OBSERVATIONS_HH_
#define AMANZI_UNSTRUCTURED_OBSERVATIONS_HH_

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Array.hpp"

#include "ObservationData.hh"
#include "Observable.hh"

#include "TimeStepManager.hh"

namespace Amanzi {

class UnstructuredObservations {

 public:
  UnstructuredObservations(Teuchos::ParameterList observations_plist,
                           const Teuchos::RCP<ObservationData>& observation_data,
                           Epetra_MpiComm* comm);

  bool DumpRequested(int cycle, double time) const;
  void MakeObservations(const State& state);
  void RegisterWithTimeStepManager(const Teuchos::Ptr<TimeStepManager>& tsm);

 private:

  typedef std::map<std::string, Teuchos::RCP<Observable> > ObservableMap;

  Teuchos::RCP<Amanzi::ObservationData> observation_data_;
  ObservableMap observations_;

};

}


#endif
