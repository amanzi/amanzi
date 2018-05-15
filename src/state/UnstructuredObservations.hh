/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/*
  ATS is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Markus Berndt
           Ethan Coon (ecoon@lanl.gov)

*/

//!  Collection of Observations on an unstructured mesh.


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
  UnstructuredObservations(Teuchos::ParameterList& observations_plist,
                           const Teuchos::RCP<ObservationData>& observation_data,
                           Epetra_MpiComm* comm);

  bool DumpRequested(int cycle, double time) const;
  void MakeObservations(const State& state);
  void RegisterWithTimeStepManager(const Teuchos::Ptr<TimeStepManager>& tsm);
  void Flush();

 private:

  typedef std::map<std::string, Teuchos::RCP<Observable> > ObservableMap;

  Teuchos::RCP<Amanzi::ObservationData> observation_data_;
  ObservableMap observations_;

};

}


#endif
