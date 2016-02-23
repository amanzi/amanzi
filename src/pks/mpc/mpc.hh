/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Interface for the Base MPC class.  A multi process coordinator is a PK
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

#ifndef PKS_MPC_MPC_HH_
#define PKS_MPC_MPC_HH_

#include <vector>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_MpiComm.h"

#include "State.hh"
#include "TreeVector.hh"

#include "pk_default_base.hh"
#include "pk_factory.hh"

namespace Amanzi {

template <class PK_t>
class MPC : virtual public PKDefaultBase {

public:
  MPC(const Teuchos::RCP<Teuchos::ParameterList>& plist,
      Teuchos::ParameterList& FElist,
      const Teuchos::RCP<TreeVector>& soln);

  // Virtual destructor
  virtual ~MPC() {}

  // PK methods
  // -- setup
  virtual void setup(const Teuchos::Ptr<State>& S);

  // -- calls all sub-PK initialize() methods
  virtual void initialize(const Teuchos::Ptr<State>& S);

  // transfer operators
  virtual void state_to_solution(const Teuchos::RCP<State>& S,
          TreeVector& soln);
  virtual void solution_to_state(TreeVector& soln,
          const Teuchos::RCP<State>& S);

  // -- loops over sub-PKs
  virtual void commit_state(double dt, const Teuchos::RCP<State>& S);
  virtual void calculate_diagnostics(const Teuchos::RCP<State>& S);

  // set States
  virtual void set_states(const Teuchos::RCP<const State>& S,
                          const Teuchos::RCP<State>& S_inter,
                          const Teuchos::RCP<State>& S_next);

 protected: // data

  typedef std::vector<Teuchos::RCP<PK_t> > SubPKList;

  SubPKList sub_pks_;

};


// -----------------------------------------------------------------------------
// Setup of PK hierarchy from PList
// -----------------------------------------------------------------------------
template <class PK_t>
MPC<PK_t>::MPC(const Teuchos::RCP<Teuchos::ParameterList>& plist,
               Teuchos::ParameterList& FElist,
               const Teuchos::RCP<TreeVector>& soln) :
    PKDefaultBase(plist, FElist, soln) {

  // loop over sub-PKs in the PK sublist, constructing the hierarchy recursively
  Teuchos::RCP<Teuchos::ParameterList> pks_list = Teuchos::sublist(plist_, "PKs");
  PKFactory pk_factory;

  if (plist_->isParameter("PKs order")) {
    // ordered
    Teuchos::Array<std::string> pk_order = plist_->get< Teuchos::Array<std::string> >("PKs order");
    int npks = pk_order.size();

    for (int i=0; i!=npks; ++i) {
      // create the solution vector
      Teuchos::RCP<TreeVector> pk_soln = Teuchos::rcp(new TreeVector());
      solution_->PushBack(pk_soln);

      // create the PK
      std::string name_i = pk_order[i];
      Teuchos::RCP<Teuchos::ParameterList> pk_list = Teuchos::sublist(pks_list,name_i);
      pk_list->set("PK name", name_i);
      Teuchos::RCP<PK_ATS> pk_notype = pk_factory.CreatePK(pk_list, FElist, pk_soln);
      Teuchos::RCP<PK_t> pk = Teuchos::rcp_dynamic_cast<PK_t>(pk_notype);
      sub_pks_.push_back(pk);
    }

  } else {
    // no order, just loop over all sublists
    for (Teuchos::ParameterList::ConstIterator i = pks_list->begin();
         i != pks_list->end(); ++i) {

      const std::string &name_i  = pks_list->name(i);
      const Teuchos::ParameterEntry  &entry_i = pks_list->entry(i);
      if (entry_i.isList()) {
        // create the solution vector
        Teuchos::RCP<TreeVector> pk_soln = Teuchos::rcp(new TreeVector());
        solution_->PushBack(pk_soln);

        // create the PK
        Teuchos::RCP<Teuchos::ParameterList> pk_list = Teuchos::sublist(pks_list,name_i);
        pk_list->set("PK name", name_i);
        Teuchos::RCP<PK_ATS> pk_notype = pk_factory.CreatePK(pk_list, FElist, pk_soln);
        Teuchos::RCP<PK_t> pk = Teuchos::rcp_dynamic_cast<PK_t>(pk_notype);
        sub_pks_.push_back(pk);
      }
    }
  }
}


// -----------------------------------------------------------------------------
// Setup of PK hierarchy from PList
// -----------------------------------------------------------------------------
template <class PK_t>
void MPC<PK_t>::setup(const Teuchos::Ptr<State>& S) {
  PKDefaultBase::setup(S);
  for (typename SubPKList::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {
    (*pk)->setup(S);
  }
}


// -----------------------------------------------------------------------------
// loop over sub-PKs, calling their initialization methods
// -----------------------------------------------------------------------------
template <class PK_t>
void MPC<PK_t>::initialize(const Teuchos::Ptr<State>& S) {
  for (typename SubPKList::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {
    (*pk)->initialize(S);
  }
};


// -----------------------------------------------------------------------------
// loop over sub-PKs, calling their state_to_solution method
// -----------------------------------------------------------------------------
template <class PK_t>
void MPC<PK_t>::state_to_solution(const Teuchos::RCP<State>& S,
                            TreeVector& soln) {
  for (unsigned int i=0; i!=sub_pks_.size(); ++i) {
    Teuchos::RCP<TreeVector> pk_soln = soln.SubVector(i);
    if (pk_soln == Teuchos::null) {
      Errors::Message message("MPC: vector structure does not match PK structure");
      Exceptions::amanzi_throw(message);
    }
    sub_pks_[i]->state_to_solution(S, *pk_soln);
  }
};


// -----------------------------------------------------------------------------
// loop over sub-PKs, calling their solution_to_state method
// -----------------------------------------------------------------------------
template <class PK_t>
void MPC<PK_t>::solution_to_state(TreeVector& soln,
                            const Teuchos::RCP<State>& S) {
  for (unsigned int i=0; i!=sub_pks_.size(); ++i) {
    Teuchos::RCP<TreeVector> pk_soln = soln.SubVector(i);
    if (pk_soln == Teuchos::null) {
      Errors::Message message("MPC: vector structure does not match PK structure");
      Exceptions::amanzi_throw(message);
    }
    sub_pks_[i]->solution_to_state(*pk_soln, S);
  }
};


// -----------------------------------------------------------------------------
// loop over sub-PKs, calling their commit_state method
// -----------------------------------------------------------------------------
template <class PK_t>
void MPC<PK_t>::commit_state(double dt, const Teuchos::RCP<State>& S) {
  for (typename SubPKList::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {
    (*pk)->commit_state(dt, S);
  }
};


// -----------------------------------------------------------------------------
// loop over sub-PKs, calling their calculate_diagnostics method
// -----------------------------------------------------------------------------
template <class PK_t>
void MPC<PK_t>::calculate_diagnostics(const Teuchos::RCP<State>& S) {
  for (typename SubPKList::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {
    (*pk)->calculate_diagnostics(S);
  }
};


// -----------------------------------------------------------------------------
// loop over sub-PKs, calling their set_states() methods
// -----------------------------------------------------------------------------
template <class PK_t>
void MPC<PK_t>::set_states(const Teuchos::RCP<const State>& S,
                     const Teuchos::RCP<State>& S_inter,
                     const Teuchos::RCP<State>& S_next) {
  PKDefaultBase::set_states(S, S_inter, S_next);

  // do the loop
  for (typename SubPKList::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {
    (*pk)->set_states(S, S_inter, S_next);
  }
};

} // close namespace Amanzi

#endif
