/* -------------------------------------------------------------------------
Amanzi

Licenses: see $ATS_DIR/COPYRIGHT, $ASCEM_DIR/COPYRIGHT
Author: Ethan Coon

Interface for the Base MPC class.  A multi process coordinator (MPC)
is a PK which coordinates other PKs.  Each of these coordinated PKs
may be MPCs themselves, or physical PKs.  Note this does NOT provide a
full implementation of PK -- it does not supply the Advance() method.
Therefore this class cannot be instantiated, but must be inherited by
derived classes which finish supplying the functionality.  Instead,
this provides the data structures and methods (which may be overridden
by derived classes) for managing multiple PKs.

Most of these methods simply loop through the coordinated PKs, calling their
respective methods.
------------------------------------------------------------------------- */

#ifndef ARCOS_MPC_HH_
#define ARCOS_MPC_HH_

#include <vector>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_MpiComm.h"

#include "State.hh"
#include "TreeVector.hh"

#include "PK.hh"
#include "PK_Factory.hh"

namespace Amanzi {

template <class PK_t>
class MPCTmp : public PK {

public:
  MPCTmp(const Teuchos::RCP<Teuchos::ParameterList>& plist,
      Teuchos::ParameterList& FElist,
      const Teuchos::RCP<TreeVector>& soln);

  // PK methods
  // -- sets up sub-PKs
  virtual void Setup(const Teuchos::Ptr<State>& S);

  // -- calls all sub-PK initialize() methods
  virtual void Initialize(const Teuchos::Ptr<State>& S);

  // -- loops over sub-PKs
  virtual void CommitState(double dt, const Teuchos::Ptr<State>& S);
  virtual void CalculateDiagnostics(const Teuchos::Ptr<State>& S);

  // -- identifier accessor
  std::string name() const { return name_; }

  // set States
  virtual void SetState(const Teuchos::RCP<State>& S);

 protected:
  // identifier
  std::string name_;

  // list of the PKs coupled by this MPC
  typedef std::vector<Teuchos::RCP<PK_t> > SubPKList;
  SubPKList sub_pks_;

  // single solution vector for the coupled problem
  Teuchos::RCP<TreeVector> solution_;

  // states
  Teuchos::RCP<State> S_;
  
};


// -----------------------------------------------------------------------------
// Setup of PK hierarchy from PList
// -----------------------------------------------------------------------------
template <class PK_t>
MPCTmp<PK_t>::MPCTmp(const Teuchos::RCP<Teuchos::ParameterList>& plist,
               Teuchos::ParameterList& FElist,
               const Teuchos::RCP<TreeVector>& soln) {
  // name the PK
  name_ = plist->name();

  // set the solution vector
  solution_ = soln;

  // loop over sub-PKs in the PK sublist, constructing the hierarchy recursively
  Teuchos::RCP<Teuchos::ParameterList> pks_list = Teuchos::sublist(plist, "PKs");
  PKFactory pk_factory;

  if (plist->isParameter("PKs order")) {
    // ordered
    Teuchos::Array<std::string> pk_order = plist->get< Teuchos::Array<std::string> >("PKs order");
    int npks = pk_order.size();

    for (int i=0; i!=npks; ++i) {
      // create the solution vector
      Teuchos::RCP<TreeVector> pk_soln = Teuchos::rcp(new TreeVector());
      solution_->PushBack(pk_soln);

      // create the PK
      std::string name_i = pk_order[i];
      Teuchos::RCP<Teuchos::ParameterList> pk_list = Teuchos::sublist(pks_list,name_i);
      pk_list->set("PK name", name_i);
      Teuchos::RCP<PK> pk_notype = pk_factory.CreatePK(pk_list, FElist, pk_soln);
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
        Teuchos::RCP<PK> pk_notype = pk_factory.CreatePK(pk_list, FElist, pk_soln);
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
void MPCTmp<PK_t>::Setup(const Teuchos::Ptr<State>& S) {
  for (typename SubPKList::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {
    (*pk)->Setup(S);
  }
}


// -----------------------------------------------------------------------------
// Loop over sub-PKs, calling their initialization methods
// -----------------------------------------------------------------------------
template <class PK_t>
void MPCTmp<PK_t>::Initialize(const Teuchos::Ptr<State>& S) {
  for (typename SubPKList::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {
    (*pk)->Initialize(S);
  }
};


// -----------------------------------------------------------------------------
// loop over sub-PKs, calling their commit state method
// -----------------------------------------------------------------------------
template <class PK_t>
void MPCTmp<PK_t>::CommitState(double dt, const Teuchos::Ptr<State>& S) {
  for (typename SubPKList::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {
    (*pk)->CommitState(dt, S);
  }
};


// -----------------------------------------------------------------------------
// loop over sub-PKs, calling their CalculateDiagnostics method
// -----------------------------------------------------------------------------
template <class PK_t>
void MPCTmp<PK_t>::CalculateDiagnostics(const Teuchos::Ptr<State>& S) {
  for (typename SubPKList::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {
    (*pk)->CalculateDiagnostics(S);
  }
};


// -----------------------------------------------------------------------------
// loop over sub-PKs, calling their set_states() methods
// -----------------------------------------------------------------------------
template <class PK_t>
void MPCTmp<PK_t>::SetState(const Teuchos::RCP<State>& S) {

  S_ = S;

  // do the loop
  for (typename SubPKList::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {
    (*pk)->SetState(S);
  }
};

} // close namespace Amanzi

#endif
