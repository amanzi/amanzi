/*
  MPC

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Walkabout
*/

#ifndef AMANZI_WALKABOUT_OBSERVATIONS_HH_
#define AMANZI_WALKABOUT_OBSERVATIONS_HH_

#include <map>

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_Array.hpp"

#include "State.hh"
#include "checkpoint.hh"

#include "TimeStepManager.hh"
#include "VerboseObject.hh"
#include "State.hh"


namespace Amanzi {

class Walkabout_observations : public Checkpoint {

public:

  Walkabout_observations(Teuchos::ParameterList& plist, Epetra_MpiComm *comm): Checkpoint (plist, comm) {};
  Walkabout_observations() : Checkpoint() {};

  void WriteWalkabout(Teuchos::RCP<State>& S);
  void CalculateDarcyVelocity(Teuchos::RCP<State>& S,
                              std::vector<AmanziGeometry::Point>& xyz, 
                              std::vector<AmanziGeometry::Point>& velocity);

  void CalculatePoreVelocity(Teuchos::RCP<State>& S,
                             std::vector<AmanziGeometry::Point>& xyz, 
                             std::vector<AmanziGeometry::Point>& velocity,
                             std::vector<double>& porosity, std::vector<double>& saturation,
                             std::vector<double>& pressure, std::vector<double>& isotherm_kd);

};

}  // namespace Amanzi

#endif
