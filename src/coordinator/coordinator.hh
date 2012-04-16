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

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Epetra_MpiComm.h"

#include "tree_vector.hh"
#include "state.hh"
#include "PK.hh"
#include "pk_factory.hh"
#include "ObservationData.H"
#include "unstructured_observations.hh"
// #include "Vis.hpp"

namespace Amanzi {

class Coordinator {

public:
  Coordinator(Teuchos::ParameterList parameter_list,
              Teuchos::RCP<Amanzi::AmanziMesh::Mesh>& mesh,
              Epetra_MpiComm* comm );
              //              Amanzi::ObservationData& output_observations);

  // PK methods
  void initialize();
  void cycle_driver();

private:
  void coordinator_init();
  void read_parameter_list();

  // PK container and factory
  Teuchos::RCP<PK> pk_;

  // states
  Teuchos::RCP<State> S_;
  Teuchos::RCP<State> S_next_;
  Teuchos::RCP<TreeVector> soln_;

  // misc setup information
  Teuchos::ParameterList parameter_list_;
  Teuchos::RCP<AmanziMesh::Mesh> mesh_;
  Teuchos::ParameterList coordinator_plist_;

  double t0_, t1_;
  double max_dt_, min_dt_;
  int end_cycle_;

  // Epetra communicator
  Epetra_MpiComm* comm_;

  // observations
  //  ObservationData& output_observations_;
  //  Teuchos::RCP<UnstructuredObservations> observations_;

  // visualization
  //Amanzi::Vis *visualization_;
};

} // close namespace Amanzi

#endif
