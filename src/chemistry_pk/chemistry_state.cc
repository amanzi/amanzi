/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

#include "chemistry_state.hh"

#include "Epetra_SerialDenseVector.h"
#include "Epetra_MultiVector.h"
#include "Teuchos_RCPDecl.hpp"

#include "State_Old.hh"

//#include "cell_geometry.hh"

#include "beaker.hh"
#include "chemistry_exception.hh"
#include "errors.hh"
#include "exceptions.hh"


namespace amanzi {
namespace chemistry {

Chemistry_State::Chemistry_State(Teuchos::RCP<State_Old> S)
    : simulation_state_(S) {

}  // end Chemistry_State


Chemistry_State::~Chemistry_State() {
}  // end ~Chemistry_State

void Chemistry_State::AllocateAdditionalChemistryStorage(
    const Beaker::BeakerComponents& components) {
  unsigned int size = components.secondary_activity_coeff.size();
  if (size > 0) {
    simulation_state_->CreateStorageSecondaryActivityCoeff(size);
  }
}  // end AllocateAdditionalChemistryStorage()

}  // namespace chemistry
}  // namespace amanzi
