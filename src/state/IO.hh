/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  State

  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  Utilities for I/O of State data.
*/

#ifndef STATE_IO_HH_
#define STATE_IO_HH_

#include <string>

// Amanzi::State
#include "Checkpoint.hh"
#include "ObservationData.hh"
#include "State.hh"
#include "Visualization.hh"

namespace Amanzi {

// Visualization
void
WriteVis(Visualization& vis, const State& S);

// Checkpointing
void
ReadCheckpoint(const Comm_ptr_type& comm, State& S, const std::string& filename);

double
ReadCheckpointInitialTime(const Comm_ptr_type& comm, std::string filename);

int
ReadCheckpointPosition(const Comm_ptr_type& comm, std::string filename);

void
ReadCheckpointObservations(const Comm_ptr_type& comm,
                           std::string filename,
                           Amanzi::ObservationData& obs_data);

void
DeformCheckpointMesh(State& S, Key domain);

// Reading from files
void
ReadVariableFromExodusII(Teuchos::ParameterList& plist, CompositeVector& var);

// Statistics
void
WriteStateStatistics(const State& S,
                     const VerboseObject& vo,
                     const Teuchos::EVerbosityLevel vl = Teuchos::VERB_HIGH);
void
WriteStateStatistics(const State& S);

} // namespace Amanzi

#endif
