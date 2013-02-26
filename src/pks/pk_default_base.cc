/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Default base with a few methods implemented in standard ways.
------------------------------------------------------------------------- */

#include "global_verbosity.hh"
#include "pk_default_base.hh"

namespace Amanzi {


void PKDefaultBase::setup(const Teuchos::Ptr<State>& S) {
  name_ = plist_.get<std::string>("PK name");

  // set up the VerboseObject
  setLinePrefix(name_);

  setDefaultVerbLevel(ATS::VerbosityLevel::level_);
  Teuchos::readVerboseObjectSublist(&plist_,this);

  // get the fancy output ??
  verbosity_ = getVerbLevel();
  out_ = getOStream();

  // cells to debug
  c0_ = plist_.get<int>("debug cell 0",0);
  c1_ = plist_.get<int>("debug cell 0",1);
}


void PKDefaultBase::set_states(const Teuchos::RCP<const State>& S,
        const Teuchos::RCP<State>& S_inter,
        const Teuchos::RCP<State>& S_next) {
  S_ = S;
  S_inter_ = S_inter;
  S_next_ = S_next;
}


void PKDefaultBase::solution_to_state(const Teuchos::RCP<const TreeVector>& soln,
        const Teuchos::RCP<State>& S) {
  Teuchos::RCP<TreeVector> nc_soln = Teuchos::rcp_const_cast<TreeVector>(soln);
  solution_to_state(nc_soln, S);
}

} // namespace
