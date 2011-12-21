/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Implementation for the Base MPC class.  A multi process coordinator is a PK
(process kernel) which coordinates several PKs.  Each of these coordinated PKs
may be MPCs themselves, or physical PKs.  Note this does NOT provide a full
implementation of PK -- it does not supply the advance() method.  Therefore
this class cannot be instantiated, but must be inherited by derived classes
which finish supplying the functionality.  Instead, this provides the data
structures and methods (which may be overridden by derived classes) for
managing multiple PKs.

Most of these methods simply loop through the coordinated PKs, calling their
respective methods.
------------------------------------------------------------------------- */

#include <string>

#include "errors.hh"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Epetra_Comm.h"
#include "Epetra_MpiComm.h"

// TODO: We are using depreciated parts of boost::filesystem
#define BOOST_FILESYSTEM_VERSION 2
#include "boost/filesystem/operations.hpp"
#include "boost/filesystem/path.hpp"

#include "dbc.hh"
#include "MPC.hh"
#include "State.hh"


namespace Amanzi {

MPC::MPC(Teuchos::ParameterList &mpc_plist,
         Teuchos::RCP<State>& S, Teuchos::RCP<TreeVector>& soln) :
    mpc_plist_(mpc_plist) {

  // do stuff for the VerboseObject
  // -- set the line prefix
  this->setLinePrefix("Amanzi::MPC         ");
  // -- make sure that the line prefix is printed
  this->getOStream()->setShowLinePrefix(true);
  // -- Read the sublist for verbosity settings.
  Teuchos::readVerboseObjectSublist(&mpc_plist_,this);

  using Teuchos::OSTab;
  Teuchos::EVerbosityLevel verbLevel = this->getVerbLevel();
  Teuchos::RCP<Teuchos::FancyOStream> out = this->getOStream();
  OSTab tab = this->getOSTab(); // This sets the line prefix and adds one tab

  // loop over sub-PKs in the PK sublist, constructing the hierarchy recursively
  Teuchos::ParameterList pks_list = mpc_plist.sublist("PKs");
  for (Teuchos::ParameterList::ConstIterator i = pks_list.begin();
       i != pks_list.end(); ++i) {

    const std::string &name_i  = pks_list.name(i);
    const Teuchos::ParameterEntry  &entry_i = pks_list.entry(i);

    if (entry_i.isList()) {
      if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_MEDIUM,true)) {
        *out << "MPC creating TreeVec with name: " << name_i << std::endl;
      }
      Teuchos::RCP<TreeVector> pk_soln = Teuchos::rcp(new TreeVector(name_i));
      soln->PushBack(pk_soln);

      if(out.get() && includesVerbLevel(verbLevel,Teuchos::VERB_MEDIUM,true)) {
        *out << "MPC creating PK with name: " << name_i << std::endl;
      }
      Teuchos::RCP<PK> pk = pk_factory_.create_pk(pks_list.sublist(name_i), S, pk_soln);
      pk->set_name(name_i);
      sub_pks_.push_back(pk);
    }
  }
};

// loop over sub-PKs, calling their initialization methods
void MPC::initialize(Teuchos::RCP<State>& S, Teuchos::RCP<TreeVector>& soln) {
  for (std::vector< Teuchos::RCP<PK> >::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {
    Teuchos::RCP<TreeVector> pk_soln;
    int subvec_not_found = soln->SubVector((*pk)->name(), pk_soln);
    if (subvec_not_found) {
      Errors::Message message("MPC: vector structure does not match PK structure");
      Exceptions::amanzi_throw(message);
    }
    (*pk)->initialize(S, pk_soln);
  }
};

// loop over sub-PKs, calling their state_to_solution method
void MPC::state_to_solution(const State& S, Teuchos::RCP<TreeVector>& soln) {
  for (std::vector< Teuchos::RCP<PK> >::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {
    Teuchos::RCP<TreeVector> pk_soln;
    int subvec_not_found = soln->SubVector((*pk)->name(), pk_soln);
    if (subvec_not_found) {
      Errors::Message message("MPC: vector structure does not match PK structure");
      Exceptions::amanzi_throw(message);
    }
    (*pk)->state_to_solution(S, pk_soln);
  }
};

// loop over sub-PKs, calling their state_to_solution method
void MPC::state_to_solution(const State& S, Teuchos::RCP<TreeVector>& soln,
                            Teuchos::RCP<TreeVector>& soln_dot) {
  for (std::vector< Teuchos::RCP<PK> >::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {
    Teuchos::RCP<TreeVector> pk_soln;
    int subvec_not_found = soln->SubVector((*pk)->name(), pk_soln);
    if (subvec_not_found) {
      Errors::Message message("MPC: vector structure does not match PK structure");
      Exceptions::amanzi_throw(message);
    }

    Teuchos::RCP<TreeVector> pk_soln_dot;
    subvec_not_found = soln_dot->SubVector((*pk)->name(), pk_soln_dot);
    if (subvec_not_found) {
      Errors::Message message("MPC: vector structure does not match PK structure");
      Exceptions::amanzi_throw(message);
    }
    (*pk)->state_to_solution(S, pk_soln, pk_soln_dot);
  }
};

// loop over sub-PKs, calling their solution_to_state method
void MPC::solution_to_state(const TreeVector& soln, Teuchos::RCP<State>& S) {
  for (std::vector< Teuchos::RCP<PK> >::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {
    Teuchos::RCP<const TreeVector> pk_soln;
    int subvec_not_found = soln.SubVector((*pk)->name(), pk_soln);
    if (subvec_not_found) {
      Errors::Message message("MPC: vector structure does not match PK structure");
      Exceptions::amanzi_throw(message);
    }
    (*pk)->solution_to_state(*pk_soln, S);
  }
};

// loop over sub-PKs, calling their solution_to_state method
void MPC::solution_to_state(const TreeVector& soln, const TreeVector& soln_dot,
                            Teuchos::RCP<State>& S) {
  for (std::vector< Teuchos::RCP<PK> >::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {
    Teuchos::RCP<const TreeVector> pk_soln;
    int subvec_not_found = soln.SubVector((*pk)->name(), pk_soln);
    if (subvec_not_found) {
      Errors::Message message("MPC: vector structure does not match PK structure");
      Exceptions::amanzi_throw(message);
    }

    Teuchos::RCP<const TreeVector> pk_soln_dot;
    subvec_not_found = soln_dot.SubVector((*pk)->name(), pk_soln_dot);
    if (subvec_not_found) {
      Errors::Message message("MPC: vector structure does not match PK structure");
      Exceptions::amanzi_throw(message);
    }
    (*pk)->solution_to_state(*pk_soln, *pk_soln_dot, S);
  }
};

// min(pk.get_dt()) for each sub-PK
double MPC::get_dt() {
  double dt = 1.e99;
  for (std::vector< Teuchos::RCP<PK> >::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {
    dt = std::min(dt, (*pk)->get_dt());
  }
  return dt;
}

// loop over sub-PKs, calling their commit_state method
void MPC::commit_state(double dt, Teuchos::RCP<State>& S) {
  for (std::vector< Teuchos::RCP<PK> >::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {
    (*pk)->commit_state(dt, S);
  }
};

// loop over sub-PKs, calling their calculate_diagnostics method
void MPC::calculate_diagnostics(Teuchos::RCP<State>& S) {
  for (std::vector< Teuchos::RCP<PK> >::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {
    (*pk)->calculate_diagnostics(S);
  }
};

// loop over sub-PKs, calling their set_states() methods
void MPC::set_states(Teuchos::RCP<const State>& S, Teuchos::RCP<State>& S_next) {
  // first set my state
  PK::set_states(S, S_next);

  // do the loop
  for (std::vector< Teuchos::RCP<PK> >::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {
    (*pk)->set_states(S, S_next);
  }
};
} // namespace
