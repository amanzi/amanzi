/*
This is the flow component of the Amanzi code. 
License: BSD
Authors: Neil Carlson (version 1) 
         Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
         Ethan Coon (ATS version) (ecoon@lanl.gov)
*/

#include "Teuchos_ParameterList.hpp"
#include "gmv_mesh.hh"

#include "flow.hh"

namespace Amanzi {
namespace Flow {

// Pointer copy of state to solution
void Flow::state_to_solution(const Teuchos::RCP<State>& S,
        const Teuchos::RCP<TreeVector>& solution) {
  solution->set_data(S->GetFieldData("pressure", "flow"));
};

// Pointer copy concentration fields from solution vector into state.  Used within
// compute_f() of strong couplers to set the current iterate in the state for
// use with other PKs.
void Flow::solution_to_state(const Teuchos::RCP<TreeVector>& solution,
        const Teuchos::RCP<State>& S) {
  S->SetData("pressure", "flow", soln->data());
};

} // namespace
} // namespace
