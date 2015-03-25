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

#ifndef _COORDINATOR_HH_
#define _COORDINATOR_HH_

#include "Teuchos_Time.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_MpiComm.h"

#include "VerboseObject.hh"

namespace Amanzi {

class TimeStepManager;
class Visualization;
class Checkpoint;
class State;
class TreeVector;
class PK;
class UnstructuredObservations;

class Coordinator {

public:
  Coordinator(Teuchos::ParameterList& parameter_list,
              Teuchos::RCP<State>& S,
              Epetra_MpiComm* comm );
              //              Amanzi::ObservationData& output_observations);

  // PK methods
  void setup();
  void initialize();
  void finalize();
  void report_memory();
  bool advance(double dt);
  void visualize(bool force=false);
  void checkpoint(double dt, bool force=false);
  double get_dt(bool after_fail=false);
  Teuchos::RCP<State> get_next_state() { return S_next_; }

  // one stop shopping
  void cycle_driver();

private:
  void coordinator_init();
  void read_parameter_list();

  // PK container and factory
  Teuchos::RCP<PK> pk_;

  // states
  Teuchos::RCP<State> S_;
  Teuchos::RCP<State> S_inter_;
  Teuchos::RCP<State> S_next_;
  Teuchos::RCP<TreeVector> soln_;

  // time step manager
  Teuchos::RCP<TimeStepManager> tsm_;

  // misc setup information
  Teuchos::RCP<Teuchos::ParameterList> parameter_list_;
  Teuchos::RCP<Teuchos::ParameterList> coordinator_list_;

  double t0_, t1_;
  double max_dt_, min_dt_;
  int cycle0_, cycle1_;

  // Epetra communicator
  Epetra_MpiComm* comm_;

  // observations
  //  ObservationData& output_observations_;
  //  Teuchos::RCP<UnstructuredObservations> observations_;

  // vis and checkpointing
  std::vector<Teuchos::RCP<Visualization> > visualization_;
  std::vector<Teuchos::RCP<Visualization> > failed_visualization_;
  Teuchos::RCP<Checkpoint> checkpoint_;
  bool restart_;
  std::string restart_filename_;

  // observations
  Teuchos::RCP<UnstructuredObservations> observations_;

  // timers
  Teuchos::RCP<Teuchos::Time> setup_timer_;
  Teuchos::RCP<Teuchos::Time> cycle_timer_;
  Teuchos::RCP<Teuchos::Time> timer_;
  double duration_;
  
  // fancy OS
  Teuchos::RCP<VerboseObject> vo_;
};

} // close namespace Amanzi

#endif
