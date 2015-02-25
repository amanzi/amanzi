/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Interface for the Coordinator.  Coordinator is basically just a class to hold
the cycle driver, which runs the overall, top level timestep loop.  It
instantiates states, ensures they are initialized, and runs the timestep loop
including Vis and restart/checkpoint dumps.  It contains one and only one PK
-- most likely this PK is an MPC of some type -- to do the actual work.
------------------------------------------------------------------------- */

#ifndef AMANZI_COORDINATOR_HH_
#define AMANZI_COORDINATOR_HH_

#include "Teuchos_Time.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_MpiComm.h"
#include "TI_Specs.hh"
#include "ObservationData.hh"
#include "Unstructured_observations.hh"

#include "VerboseObject.hh"

namespace Amanzi {

class TimeStepManager;
class Visualization;
class Checkpoint;
class State;
class TreeVector;
class PK;
class UnstructuredObservations;

class CycleDriver {

public:
  CycleDriver(Teuchos::RCP<Teuchos::ParameterList> glist_,
              Teuchos::RCP<State>& S,
              Epetra_MpiComm* comm,
              Amanzi::ObservationData& output_observations);


  // PK methods
  void setup();
  void initialize();
  void init_pk(int);
  void reset_pk();
  void finalize();
  void report_memory();
  bool advance(double dt);
  void visualize(bool force=false);
  void observations(bool force=false);
  void checkpoint(double dt, bool force=false);
  double get_dt();
  void set_dt(double dt);
  void reset_driver(int time_period_id);
  // one stop shopping
  void go();

private:
  void coordinator_init_();
  void read_parameter_list_();

  // PK container and factory
  Teuchos::RCP<PK> pk_;

  // states
  Teuchos::RCP<State> S_, S_old_;
  Teuchos::RCP<TreeVector> soln_;

  // time step manager
  Teuchos::Ptr<TimeStepManager> tsm_;

  // misc setup information
  Teuchos::RCP<Teuchos::ParameterList> parameter_list_;
  Teuchos::RCP<Teuchos::ParameterList> coordinator_list_;

  double t0_, t1_;
  std::vector<double> t_, tp_start_, tp_end_, tp_dt_, tp_max_cycle_;  
  double max_dt_, min_dt_;
  int cycle0_, cycle1_;
  int num_time_periods_;
  int time_period_id_;

  // Epetra communicator
  Epetra_MpiComm* comm_;

  // observations
  Amanzi::ObservationData&  output_observations_;
  Teuchos::RCP<Amanzi::Unstructured_observations> observations_;

  // time interval
  std::vector<Amanzi::Flow::TI_Specs> ti_specs_;


  // vis and checkpointing
  std::vector<Teuchos::RCP<Visualization> > visualization_;
  std::vector<Teuchos::RCP<Visualization> > failed_visualization_;
  Teuchos::RCP<Checkpoint> checkpoint_;
  bool restart_requested_;
  std::string restart_filename_;

  // time period control
  std::vector<std::pair<double,double> > reset_info_;

  // //  checkpoint/restart 
  // Teuchos::RCP<Amanzi::Checkpoint> restart_;
 
  // walkabout
  Teuchos::RCP<Amanzi::Checkpoint> walkabout;  

  // fancy OS
  Teuchos::RCP<VerboseObject> vo_;
};

} // close namespace Amanzi

#endif
