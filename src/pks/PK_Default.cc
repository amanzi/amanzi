/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon, Daniil Svyatsky

Default base with a few methods implemented in standard ways.
------------------------------------------------------------------------- */

#include "PK_Default.hh"

namespace Amanzi {


void PK_Default::Setup(const Teuchos::Ptr<State>& S) {

  // // THIS MAY BE CALLED MORE THAN ONCE!
  // name_ = plist_->get<std::string>("PK name");
  // name_ = pk_tree_.name()

  // set up the VerboseObject
  // vo_ = Teuchos::rcp(new VerboseObject(name_, *plist_));
}


void PK_Default::set_states(const Teuchos::RCP<const State>& S,
                              const Teuchos::RCP<State>& S_inter,
                              const Teuchos::RCP<State>& S_next) {
    S_ = S;
    S_inter_ = S_inter;
    S_next_ = S_next;
}


void PK_Default::Solution_to_State(const TreeVector& soln,
        const Teuchos::RCP<State>& S) {
  TreeVector* soln_nc_ptr = const_cast<TreeVector*>(&soln);
  Solution_to_State(*soln_nc_ptr, S);
}

} // namespace
