/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

/*
  MPC

  Checkpointing Walkabout data.
*/

#ifndef AMANZI_WALKABOUT_CHECKPOINT_HH_
#define AMANZI_WALKABOUT_CHECKPOINT_HH_

#include <map>

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Array.hpp"

#include "PK.hh"
#include "State.hh"
#include "Checkpoint.hh"

#include "TimeStepManager.hh"
#include "VerboseObject.hh"


namespace Amanzi {

class WalkaboutCheckpoint : public Checkpoint {
 public:
  WalkaboutCheckpoint(Teuchos::ParameterList& plist, const State& S) : Checkpoint(plist, S){};
  WalkaboutCheckpoint() : Checkpoint(true){};

  // output of fields
  void WriteDataFile(Teuchos::RCP<State>& S, Teuchos::RCP<PK> pk);

  // recontruct vector velocity at mesh nodes
  void CalculateDarcyVelocity(Teuchos::RCP<State>& S,
                              std::vector<AmanziGeometry::Point>& xyz,
                              std::vector<AmanziGeometry::Point>& velocity) const;

  // interpolate various fileds to mesh nodes
  void CalculateData(Teuchos::RCP<State>& S,
                     std::vector<AmanziGeometry::Point>& xyz,
                     std::vector<AmanziGeometry::Point>& velocity,
                     std::vector<double>& porosity,
                     std::vector<double>& saturation,
                     std::vector<double>& pressure,
                     std::vector<double>& isotherm_kd,
                     std::vector<int>& material_ids);

 private:
  Teuchos::RCP<PK> pk_;
};

} // namespace Amanzi

#endif
