/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

#include "chemistry_state.hh"

#include "Epetra_SerialDenseVector.h"
#include "Epetra_MultiVector.h"
#include "Teuchos_RCPDecl.hpp"

#include "State.hpp"

#include "cell_geometry.hh"

#include "chemistry_exception.hh"
#include "errors.hh"
#include "exceptions.hh"


namespace amanzi {
namespace chemistry {

Chemistry_State::Chemistry_State(Teuchos::RCP<State> S)
    : simulation_state_(S) {

}  // end Chemistry_State


Chemistry_State::~Chemistry_State() {
}  // end ~Chemistry_State


}  // namespace chemistry
}  // namespace amanzi
