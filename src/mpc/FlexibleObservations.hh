/*
  Multi-Process Coordinator

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Flexible unstructured observations.
*/

#ifndef AMANZI_FLEXIBLE_OBSERVATIONS_HH_
#define AMANZI_FLEXIBLE_OBSERVATIONS_HH_

#include <map>

// TPLs
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Array.hpp"

// Amanzi
#include "IOEvent.hh"
#include "ObservationData.hh"
#include "State.hh"
#include "TimeStepManager.hh"
#include "Units.hh"
#include "VerboseObject.hh"
#include "ObservableAmanzi.hh"
#include "Key.hh"

namespace Amanzi {

class FlexibleObservations {
 public:
  FlexibleObservations(Teuchos::RCP<Teuchos::ParameterList> coordinator_list,
                       Teuchos::RCP<Teuchos::ParameterList> obs_list,
                       Teuchos::RCP<Teuchos::ParameterList> units_list,
                       Amanzi::ObservationData& observation_data,
                       Teuchos::RCP<const State> S);

  ~FlexibleObservations()
  {
    if (vo_ != NULL) delete vo_;
  }

  void RegisterComponentNames(std::vector<std::string>& comp_names,
                              std::vector<double>& comp_mol_masses,
                              int num_liquid)
  {
    comp_names_ = comp_names;
    comp_mol_masses_ = comp_mol_masses;
    num_liquid_ = num_liquid;
  }

  int MakeObservations(State& S);
  int MakeContinuousObservations(State& S);

  bool DumpRequested(const int);
  bool DumpRequested(const double);
  bool DumpRequested(const int, const double);

  void RegisterWithTimeStepManager(const Teuchos::Ptr<TimeStepManager>& tsm);

  void FlushObservations();

 private:
  double CalculateWaterTable_(State& S, AmanziMesh::Entity_ID_List& ids);

 protected:
  VerboseObject* vo_;

 private:
  int rank_;
  Teuchos::RCP<Teuchos::ParameterList> obs_list_;
  Teuchos::RCP<Teuchos::ParameterList> coordinator_list_;
  Amanzi::ObservationData& observation_data_;
  std::map<std::string, Teuchos::RCP<Observable>> observations;

  std::vector<std::string> comp_names_;
  std::vector<double> comp_mol_masses_;
  int num_liquid_;

  Utils::Units units_;
};

} // namespace Amanzi

#endif
