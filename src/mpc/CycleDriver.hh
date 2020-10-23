/*
  Multi-Process Coordinator

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
           Daniil Svyatskiy

  Interface for the Coordinator.  Coordinator is basically just a class to hold
  the cycle driver, which runs the overall, top level timestep loop.  It
  instantiates states, ensures they are initialized, and runs the timestep loop
  including Vis and restart/checkpoint dumps.  It contains one and only one PK
  -- most likely this PK is an MPC of some type -- to do the actual work.
*/

#ifndef AMANZI_CYCLE_DRIVER_HH_
#define AMANZI_CYCLE_DRIVER_HH_

#include "Epetra_MpiComm.h"
#include "Teuchos_Time.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "FlexibleObservations.hh"
#include "ObservationData.hh"
#include "WalkaboutCheckpoint.hh"
#include "VerboseObject.hh"

namespace Amanzi {

class TimeStepManager;
class Visualization;
class Checkpoint;
class State;
class TreeVector;
class PK;

class CycleDriver {
 public:
  CycleDriver(Teuchos::RCP<Teuchos::ParameterList> glist,
              Teuchos::RCP<Amanzi::State>& S,
              const Comm_ptr_type& comm,
              Amanzi::ObservationData& observations_data);

  // PK methods
  void Setup();
  void Initialize();
  void Init_PK(int);
  void Reset_PK();
  void Finalize();
  void ReportMemory();
  double Advance(double dt);
  void Visualize(bool force = false, const std::string& tag = "");
  void Observations(bool force = false);
  void WriteCheckpoint(double dt, bool force = false);
  void WriteWalkabout(bool force);
  //  void RegisterOutput();
  double get_dt(bool after_failuer = false);
  void set_dt(double dt);
  void ResetDriver(int time_period_id);
  // one stop shopping
  Teuchos::RCP<State> Go();

  // access (for unit tests only)
  Teuchos::RCP<const Amanzi::WalkaboutCheckpoint> walkabout() const { return walkabout_; }

 private:
  void CoordinatorInit_();
  void ReadParameterList_();

 private:
  // PK container and factory
  Teuchos::RCP<PK> pk_;

  // states
  Teuchos::RCP<State> S_, S_old_;
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  Teuchos::RCP<TreeVector> soln_;

  // time step manager
  Teuchos::RCP<TimeStepManager> tsm_;

  // misc setup information
  Teuchos::RCP<Teuchos::ParameterList> glist_;
  Teuchos::RCP<Teuchos::ParameterList> coordinator_list_;

  std::vector<double> t_, tp_start_, tp_end_, tp_dt_, tp_max_cycle_, tp_max_dt_;  
  double max_dt_, min_dt_;
  int cycle0_, cycle1_;
  int num_time_periods_;
  int time_period_id_;

  // Amanzi communicator
  Comm_ptr_type comm_;

  // observations
  Amanzi::ObservationData& observations_data_;
  Teuchos::RCP<FlexibleObservations> observations_;

  // vis and checkpointing
  std::vector<Teuchos::RCP<Visualization> > visualization_;
  std::vector<Teuchos::RCP<Visualization> > failed_visualization_;
  Teuchos::RCP<Checkpoint> checkpoint_;
  bool restart_requested_;
  //  bool output_registered_;
  std::string restart_filename_;

  // time period control
  std::vector<std::pair<double,double> > reset_info_;
  std::vector<std::pair<double,double> > reset_max_;

  // //  checkpoint/restart 
  // Teuchos::RCP<Amanzi::Checkpoint> restart_;
 
  // walkabout
  Teuchos::RCP<Amanzi::WalkaboutCheckpoint> walkabout_;  

  // fancy OS
  Teuchos::RCP<VerboseObject> vo_;
};


// non-meber function
// names of all fields that go into a vis file
inline
std::set<std::string> StateVisFields(const State& S)
{
  std::set<std::string> fields;
  for (auto it = S.field_begin(); it != S.field_end(); ++it) 
    if (it->second->io_vis()) fields.insert(it->first);
  return fields;
}

}  // namespace Amanzi

#endif
