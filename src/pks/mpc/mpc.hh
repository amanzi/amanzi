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
#include "Teuchos_VerboseObject.hpp"
#include "Epetra_MpiComm.h"

#include "state.hh"
#include "tree_vector.hh"

#include "pk_default_base.hh"
#include "pk_factory.hh"

namespace Amanzi {

template <class PK_Type>
class MPC : virtual public PKDefaultBase {

public:
  MPC(Teuchos::ParameterList& mpc_plist, const Teuchos::RCP<TreeVector>& soln) :
      PKDefaultBase(mpc_plist,soln) {}

  // Virtual destructor
  virtual ~MPC() {}

  // PK methods
  // -- setup
  virtual void setup(const Teuchos::Ptr<State>& S);

  // -- calls all sub-PK initialize() methods
  virtual void initialize(const Teuchos::Ptr<State>& S);

  // transfer operators
  virtual void state_to_solution(const Teuchos::RCP<State>& S,
          const Teuchos::RCP<TreeVector>& soln);
  virtual void solution_to_state(const Teuchos::RCP<TreeVector>& soln,
          const Teuchos::RCP<State>& S);

  // -- loops over sub-PKs
  virtual void commit_state(double dt, const Teuchos::RCP<State>& S);
  virtual void calculate_diagnostics(const Teuchos::RCP<State>& S);

  // set States
  virtual void set_states(const Teuchos::RCP<const State>& S,
                          const Teuchos::RCP<State>& S_inter,
                          const Teuchos::RCP<State>& S_next);

 protected: // data

  typedef std::vector<Teuchos::RCP<PK_Type> > SubPKList;

  SubPKList sub_pks_;

};

// -----------------------------------------------------------------------------
// Setup of PK hierarchy from PList
// -----------------------------------------------------------------------------
template <class PK_Type>
void MPC<PK_Type>::setup(const Teuchos::Ptr<State>& S) {
  PKDefaultBase::setup(S);

  // loop over sub-PKs in the PK sublist, constructing the hierarchy recursively
  Teuchos::ParameterList pks_list = plist_.sublist("PKs");
  PKFactory pk_factory;

  if (plist_.isParameter("PKs order")) {
    // ordered
    Teuchos::Array<std::string> pk_order = plist_.get< Teuchos::Array<std::string> >("PKs order");
    int npks = pk_order.size();

    for (int i=0; i!=npks; ++i) {
      // create the solution vector
      std::string name_i = pk_order[i];
      Teuchos::RCP<TreeVector> pk_soln = Teuchos::rcp(new TreeVector(name_i));
      solution_->PushBack(pk_soln);

      // create the PK
      Teuchos::ParameterList pk_list = pks_list.sublist(name_i);
      pk_list.set("PK name", name_i);
      Teuchos::RCP<PK> pk_notype = pk_factory.CreatePK(pk_list, pk_soln);
      pk_notype->setup(S);

      Teuchos::RCP<PK_Type> pk = Teuchos::rcp_dynamic_cast<PK_Type>(pk_notype);
      sub_pks_.push_back(pk);
    }

  } else {
    // no order, just loop over all sublists
    for (Teuchos::ParameterList::ConstIterator i = pks_list.begin();
         i != pks_list.end(); ++i) {

      const std::string &name_i  = pks_list.name(i);
      const Teuchos::ParameterEntry  &entry_i = pks_list.entry(i);
      if (entry_i.isList()) {
        // create the solution vector
        Teuchos::RCP<TreeVector> pk_soln = Teuchos::rcp(new TreeVector(name_i));
        solution_->PushBack(pk_soln);

        // create the PK
        Teuchos::ParameterList pk_list = pks_list.sublist(name_i);
        pk_list.set("PK name", name_i);
        Teuchos::RCP<PK> pk_notype = pk_factory.CreatePK(pk_list, pk_soln);
        pk_notype->setup(S);

        Teuchos::RCP<PK_Type> pk = Teuchos::rcp_dynamic_cast<PK_Type>(pk_notype);
        sub_pks_.push_back(pk);
      }
    }
  }
}


// -----------------------------------------------------------------------------
// loop over sub-PKs, calling their initialization methods
// -----------------------------------------------------------------------------
template <class PK_Type>
void MPC<PK_Type>::initialize(const Teuchos::Ptr<State>& S) {
  for (typename SubPKList::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {
    (*pk)->initialize(S);
  }
};


// -----------------------------------------------------------------------------
// loop over sub-PKs, calling their state_to_solution method
// -----------------------------------------------------------------------------
template <class PK_Type>
void MPC<PK_Type>::state_to_solution(const Teuchos::RCP<State>& S,
                            const Teuchos::RCP<TreeVector>& soln) {
  for (typename SubPKList::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {
    Teuchos::RCP<TreeVector> pk_soln = soln->SubVector((*pk)->name());
    if (pk_soln == Teuchos::null) {
      Errors::Message message("MPC: vector structure does not match PK structure");
      Exceptions::amanzi_throw(message);
    }
    (*pk)->state_to_solution(S, pk_soln);
  }
};


// -----------------------------------------------------------------------------
// loop over sub-PKs, calling their solution_to_state method
// -----------------------------------------------------------------------------
template <class PK_Type>
void MPC<PK_Type>::solution_to_state(const Teuchos::RCP<TreeVector>& soln,
                            const Teuchos::RCP<State>& S) {
  for (typename SubPKList::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {
    Teuchos::RCP<TreeVector> pk_soln = soln->SubVector((*pk)->name());
    if (pk_soln == Teuchos::null) {
      Errors::Message message("MPC: vector structure does not match PK structure");
      Exceptions::amanzi_throw(message);
    }
    (*pk)->solution_to_state(pk_soln, S);
  }
};


// -----------------------------------------------------------------------------
// loop over sub-PKs, calling their commit_state method
// -----------------------------------------------------------------------------
template <class PK_Type>
void MPC<PK_Type>::commit_state(double dt, const Teuchos::RCP<State>& S) {
  for (typename SubPKList::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {
    (*pk)->commit_state(dt, S);
  }
};


// -----------------------------------------------------------------------------
// loop over sub-PKs, calling their calculate_diagnostics method
// -----------------------------------------------------------------------------
template <class PK_Type>
void MPC<PK_Type>::calculate_diagnostics(const Teuchos::RCP<State>& S) {
  for (typename SubPKList::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {
    (*pk)->calculate_diagnostics(S);
  }
};


// -----------------------------------------------------------------------------
// loop over sub-PKs, calling their set_states() methods
// -----------------------------------------------------------------------------
template <class PK_Type>
void MPC<PK_Type>::set_states(const Teuchos::RCP<const State>& S,
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
