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

//#include "pk_default_base.hh"
//#include "pk_factory_ats.hh"
#include "PK.hh"
#include "PK_Factory.hh"

namespace Amanzi {

template <class PK_t>
class MPC : virtual public PK {

public:
  MPC(Teuchos::ParameterList& FElist,
      const Teuchos::RCP<Teuchos::ParameterList>& plist,
      const Teuchos::RCP<State>& S,
      const Teuchos::RCP<TreeVector>& solution);

  // Virtual destructor
  virtual ~MPC() {}

  // PK methods
  // -- setup
  virtual void Setup(const Teuchos::Ptr<State>& S);

  // -- calls all sub-PK initialize() methods
  virtual void Initialize(const Teuchos::Ptr<State>& S);

  // transfer operators
  virtual void State_to_Solution(const Teuchos::RCP<State>& S,
          TreeVector& soln);
  virtual void Solution_to_State(TreeVector& soln,
          const Teuchos::RCP<State>& S);
  virtual void Solution_to_State(const TreeVector& soln,
          const Teuchos::RCP<State>& S);

  // -- loops over sub-PKs
  virtual void CommitStep(double t_old, double t_new, const Teuchos::RCP<State>& S);
  virtual void CalculateDiagnostics(const Teuchos::RCP<State>& S);
  virtual bool valid_step();
  
  // set States
  virtual void set_states(const Teuchos::RCP<const State>& S,
                          const Teuchos::RCP<State>& S_inter,
                          const Teuchos::RCP<State>& S_next);

 protected: // data

  typedef std::vector<Teuchos::RCP<PK_t> > SubPKList;
  Teuchos::RCP<Teuchos::ParameterList> my_list_;
  Teuchos::RCP<Teuchos::ParameterList> global_list_;
  SubPKList sub_pks_;

};


// -----------------------------------------------------------------------------
// Setup of PK hierarchy from PList
// -----------------------------------------------------------------------------
template <class PK_t>
MPC<PK_t>::MPC(Teuchos::ParameterList& FElist,
      const Teuchos::RCP<Teuchos::ParameterList>& plist,
      const Teuchos::RCP<State>& S,
      const Teuchos::RCP<TreeVector>& solution)
//    PKDefaultBase(plist, FElist, soln)
{

  plist_ = plist;
  solution_ = solution;

  // loop over sub-PKs in the PK sublist, constructing the hierarchy recursively
  Teuchos::RCP<Teuchos::ParameterList> pks_list = Teuchos::sublist(plist_, "PKs");;

  my_list_ = Teuchos::sublist(plist_, "PKs");
  
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
      //std::cout << *pk_list<<"\n";
      pk_list->set("PK name", name_i);
      const std::string &pk_origin = pk_list -> get<std::string>("PK origin", "ATS");
      //      std::cout<<*pk_list<<"\n";
      if (pk_origin == "ATS"){
        Teuchos::RCP<PK> pk_notype = pk_factory.CreatePK(FElist, pk_list, S, pk_soln);
        Teuchos::RCP<PK_t> pk = Teuchos::rcp_dynamic_cast<PK_t>(pk_notype);
        sub_pks_.push_back(pk);
      }
      else if (pk_origin == "Amanzi"){
        Teuchos::RCP<PK> pk_notype = pk_factory.CreatePK(*pk_list, plist_, S, pk_soln);
        Teuchos::RCP<PK_t> pk = Teuchos::rcp_dynamic_cast<PK_t>(pk_notype);
        sub_pks_.push_back(pk);
      }
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
        //std::cout << *pk_list<<"\n";
        pk_list->set("PK name", name_i);
        const std::string &pk_origin = pk_list -> get<std::string>("PK origin", "ATS");
        
        if (pk_origin == "ATS"){
          Teuchos::RCP<PK> pk_notype = pk_factory.CreatePK(FElist, pk_list, S, pk_soln);
          Teuchos::RCP<PK_t> pk = Teuchos::rcp_dynamic_cast<PK_t>(pk_notype);
          sub_pks_.push_back(pk);
        }
        else if (pk_origin == "Amanzi"){
          Teuchos::RCP<PK> pk_notype = pk_factory.CreatePK(*pk_list, plist_, S, pk_soln);
          Teuchos::RCP<PK_t> pk = Teuchos::rcp_dynamic_cast<PK_t>(pk_notype);
          sub_pks_.push_back(pk);
        }
      }
    }
  }
}


// -----------------------------------------------------------------------------
// Setup of PK hierarchy from PList
// -----------------------------------------------------------------------------
template <class PK_t>
void MPC<PK_t>::Setup(const Teuchos::Ptr<State>& S) {
  //PKDefaultBase::setup(S);
  // THIS MAY BE CALLED MORE THAN ONCE!
  name_ = plist_->get<std::string>("PK name");
  // set up the VerboseObject
  vo_ = Teuchos::rcp(new VerboseObject(name_, *plist_));

  for (typename SubPKList::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {
    (*pk)->Setup(S);
  }
}


// -----------------------------------------------------------------------------
// loop over sub-PKs, calling their initialization methods
// -----------------------------------------------------------------------------
template <class PK_t>
void MPC<PK_t>::Initialize(const Teuchos::Ptr<State>& S) {
  for (typename SubPKList::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {
    (*pk)->Initialize(S);
  }
};


// -----------------------------------------------------------------------------
// loop over sub-PKs, calling their state_to_solution method
// -----------------------------------------------------------------------------
template <class PK_t>
void MPC<PK_t>::State_to_Solution(const Teuchos::RCP<State>& S,
                            TreeVector& soln) {
  for (unsigned int i=0; i!=sub_pks_.size(); ++i) {
    Teuchos::RCP<TreeVector> pk_soln = soln.SubVector(i);
    if (pk_soln == Teuchos::null) {
      Errors::Message message("MPC: vector structure does not match PK structure");
      Exceptions::amanzi_throw(message);
    }
    sub_pks_[i]->State_to_Solution(S, *pk_soln);
  }
};


// -----------------------------------------------------------------------------
// loop over sub-PKs, calling their solution_to_state method
// -----------------------------------------------------------------------------
template <class PK_t>
void MPC<PK_t>::Solution_to_State(TreeVector& soln,
                            const Teuchos::RCP<State>& S) {
  for (unsigned int i=0; i!=sub_pks_.size(); ++i) {
    Teuchos::RCP<TreeVector> pk_soln = soln.SubVector(i);
    if (pk_soln == Teuchos::null) {
      Errors::Message message("MPC: vector structure does not match PK structure");
      Exceptions::amanzi_throw(message);
    }
    sub_pks_[i]->Solution_to_State(*pk_soln, S);
  }
};

template <class PK_t>
void MPC<PK_t>::Solution_to_State(const TreeVector& soln,
                            const Teuchos::RCP<State>& S) {
  TreeVector* soln_nc_ptr = const_cast<TreeVector*>(&soln);
  Solution_to_State(*soln_nc_ptr, S);
}


// -----------------------------------------------------------------------------
// loop over sub-PKs, calling their commit_state method
// -----------------------------------------------------------------------------
template <class PK_t>
void MPC<PK_t>::CommitStep(double t_old, double t_new, const Teuchos::RCP<State>& S) {
  for (typename SubPKList::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {
    (*pk)->CommitStep(t_old, t_new, S);
  }
};


// -----------------------------------------------------------------------------
// loop over sub-PKs, calling their calculate_diagnostics method
// -----------------------------------------------------------------------------
template <class PK_t>
void MPC<PK_t>::CalculateDiagnostics(const Teuchos::RCP<State>& S) {
  for (typename SubPKList::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {
    (*pk)->CalculateDiagnostics(S);
  }
};


// -----------------------------------------------------------------------------
// loop over sub-PKs, calling their valid_step() method
// -----------------------------------------------------------------------------
template <class PK_t>
bool MPC<PK_t>::valid_step() {
  bool valid(true);
  // for (typename SubPKList::iterator pk = sub_pks_.begin();
  //      pk != sub_pks_.end(); ++pk) {
  //   bool is_valid = (*pk)->valid_step();
  //   if (!is_valid) valid = is_valid;
  // }
  return valid;
};


// -----------------------------------------------------------------------------
// loop over sub-PKs, calling their set_states() methods
// -----------------------------------------------------------------------------
template <class PK_t>
void MPC<PK_t>::set_states(const Teuchos::RCP<const State>& S,
                     const Teuchos::RCP<State>& S_inter,
                     const Teuchos::RCP<State>& S_next) {
  //  PKDefaultBase::set_states(S, S_inter, S_next);
  S_ = S;
  S_inter_ = S_inter;
  S_next_ = S_next;

  //  do the loop
  for (typename SubPKList::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {
    (*pk)->set_states(S, S_inter, S_next);
  }
};

} // close namespace Amanzi

#endif
