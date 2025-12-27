/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
           Daniil Svyatskiy
*/

/*
  MPC PK

  Interface for the Base MPC class.  A multi process coordinator (MPC)
  is a PK which coordinates other PKs.  Each of these coordinated PKs
  may be MPCs themselves, or physical PKs.  Note this does NOT provide a
  full implementation of PK -- it does not supply the AdvanceStep() method.
  Therefore this class cannot be instantiated, but must be inherited by
  derived classes which finish supplying the functionality.  Instead,
  this provides the data structures and methods (which may be overridden
  by derived classes) for managing multiple PKs.

  Most of these methods simply loop through the coordinated PKs, calling their
  respective methods.
*/

#ifndef AMANZI_PK_MPC_HH_
#define AMANZI_PK_MPC_HH_

#include <vector>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_MpiComm.h"

#include "LScheme_Helpers.hh"
#include "State.hh"
#include "TreeVector.hh"

#include "PK.hh"
#include "PK_Factory.hh"

namespace Amanzi {

template<class PK_Base>
class PK_MPC : virtual public PK {
 public:
  PK_MPC() {};
  PK_MPC(Teuchos::ParameterList& pk_tree,
         const Teuchos::RCP<Teuchos::ParameterList>& global_list,
         const Teuchos::RCP<State>& S,
         const Teuchos::RCP<TreeVector>& soln);
  ~PK_MPC() {};

  // PK methods
  // -- sets up sub-PKs
  virtual void parseParameterList();
  virtual void Setup();
  virtual void Initialize();

  // -- sets L-scheme values for sub-PKs 
  virtual std::vector<Key> SetupLSchemeKey(Teuchos::ParameterList& plist);
  virtual void InitializeLSchemeStep();
  virtual void ComputeLSchemeStability();

  // -- set tags for integration period
  virtual void set_tags(const Tag& current, const Tag& next);

  // -- loops over sub-PKs
  virtual void CommitStep(double t_old, double t_new, const Tag& tag);
  virtual void CalculateDiagnostics(const Tag& tag);

  virtual void State_to_Solution(const Tag& tag, TreeVector& soln) {};
  virtual void Solution_to_State(const TreeVector& soln, const Tag& tag) {};

  // iterator over pks
  typename std::vector<Teuchos::RCP<PK_Base>>::iterator begin() { return sub_pks_.begin(); }
  typename std::vector<Teuchos::RCP<PK_Base>>::iterator end() { return sub_pks_.end(); }

  Teuchos::RCP<PK_Base> getSubPK(int i) { return sub_pks_[i]; }

 protected:
  Teuchos::RCP<Teuchos::ParameterList> getSubPKPlist_(int i);

  // list of the PKs coupled by this MPC
  typedef std::vector<Teuchos::RCP<PK_Base>> SubPKList;
  SubPKList sub_pks_;

  // single solution vector for this pk only
  Teuchos::RCP<TreeVector> my_solution_;

  // L-scheme data
  bool L_scheme_ = false;
  std::vector<Key> L_scheme_keys_;

  // plists
  Teuchos::RCP<Teuchos::ParameterList> global_list_;
  Teuchos::RCP<Teuchos::ParameterList> my_list_;
  Teuchos::ParameterList pk_tree_;
};


// -----------------------------------------------------------------------------
// Costructor
// -----------------------------------------------------------------------------
template<class PK_Base>
PK_MPC<PK_Base>::PK_MPC(Teuchos::ParameterList& pk_tree,
                        const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                        const Teuchos::RCP<State>& S,
                        const Teuchos::RCP<TreeVector>& soln)
  : global_list_(global_list), pk_tree_(pk_tree)
{
  S_ = S;

  // instead of calling the base class contructor, we initialize here
  solution_ = soln;

  // name the PK
  name_ = pk_tree.name();
  auto found = name_.rfind("->");
  if (found != std::string::npos) name_.erase(0, found + 2);

  // get my parameter list
  my_list_ = Teuchos::sublist(Teuchos::sublist(global_list_, "PKs"), name_);

  Teuchos::RCP<Teuchos::ParameterList> plist;
  if (global_list_->isSublist("PKs")) {
    plist = Teuchos::sublist(global_list, "PKs");
  }

  std::vector<std::string> pk_name =
    my_list_->get<Teuchos::Array<std::string>>("PKs order").toVector();

  // loop over sub-PKs in the PK sublist, constructing the hierarchy recursively
  PKFactory pk_factory;

  for (int i = 0; i < pk_name.size(); ++i) {
    const std::string& sub_name = pk_name[i];
    if (!plist->isSublist(sub_name)) {
      Errors::Message msg;
      msg << "PK tree has no sublist \"" << sub_name << "\"";
      Exceptions::amanzi_throw(msg);
    }
  }

  my_solution_ = Teuchos::rcp(new TreeVector());
  sub_pks_.clear();
  for (int i = 0; i < pk_name.size(); i++) {
    // create vector of primary unknowns
    Teuchos::RCP<TreeVector> pk_soln = Teuchos::rcp(new TreeVector());
    solution_->PushBack(pk_soln);
    my_solution_->PushBack(pk_soln);

    // create the PK
    Teuchos::RCP<PK> pk_notype = pk_factory.CreatePK(pk_name[i], pk_tree, global_list, S, pk_soln);
    Teuchos::RCP<PK_Base> pk = Teuchos::rcp_dynamic_cast<PK_Base>(pk_notype);
    sub_pks_.push_back(pk);
  }

  // share names with siblings
  for (int n = 0; n < pk_name.size(); ++n) {
    for (int i = 0; i < pk_name.size(); ++i) {
      if (i != n) {
        Teuchos::RCP<PK> pk = Teuchos::rcp_dynamic_cast<PK>(sub_pks_[i]);
        plist->sublist(pk_name[n]).set<Teuchos::RCP<PK>>("sibling " + std::to_string(i), pk);
      }
    }
  }
}


// -----------------------------------------------------------------------------
// Parse parameter list.
// -----------------------------------------------------------------------------
template<class PK_Base>
void
PK_MPC<PK_Base>::parseParameterList()
{
  for (typename SubPKList::iterator pk = sub_pks_.begin(); pk != sub_pks_.end(); ++pk) {
    (*pk)->parseParameterList();
  }
}


// -----------------------------------------------------------------------------
// Setup of PK hierarchy from PList
// -----------------------------------------------------------------------------
template<class PK_Base>
void
PK_MPC<PK_Base>::Setup()
{
  for (typename SubPKList::iterator pk = sub_pks_.begin(); pk != sub_pks_.end(); ++pk) {
    (*pk)->Setup();
  }
}


// -----------------------------------------------------------------------------
// Loop over sub-PKs, calling their initialization methods
// -----------------------------------------------------------------------------
template<class PK_Base>
void
PK_MPC<PK_Base>::Initialize()
{
  for (typename SubPKList::iterator pk = sub_pks_.begin(); pk != sub_pks_.end(); ++pk) {
    (*pk)->Initialize();
  }
}


// -----------------------------------------------------------------------------
// Propagate tags to sub-PKs
// -----------------------------------------------------------------------------
template<class PK_Base>
void
PK_MPC<PK_Base>::set_tags(const Tag& current, const Tag& next)
{
  PK::set_tags(current, next);
  for (auto& pk : sub_pks_) pk->set_tags(current, next);
}


// -----------------------------------------------------------------------------
// loop over sub-PKs, calling their commit state method
// -----------------------------------------------------------------------------
template<class PK_Base>
void
PK_MPC<PK_Base>::CommitStep(double t_old, double t_new, const Tag& tag)
{
  for (typename SubPKList::iterator pk = sub_pks_.begin(); pk != sub_pks_.end(); ++pk) {
    (*pk)->CommitStep(t_old, t_new, tag);
  }
}


// -----------------------------------------------------------------------------
// loop over sub-PKs, calling their CalculateDiagnostics method
// -----------------------------------------------------------------------------
template<class PK_Base>
void
PK_MPC<PK_Base>::CalculateDiagnostics(const Tag& tag)
{
  for (typename SubPKList::iterator pk = sub_pks_.begin(); pk != sub_pks_.end(); ++pk) {
    (*pk)->CalculateDiagnostics(tag);
  }
}


template<class PK_Base>
Teuchos::RCP<Teuchos::ParameterList>
PK_MPC<PK_Base>::getSubPKPlist_(int i)
{
  Teuchos::Array<std::string> names = my_list_->get<Teuchos::Array<std::string>>("PKs order");
  return Teuchos::sublist(Teuchos::sublist(global_list_, "PKs"), names[i]);
}


// -----------------------------------------------------------------------------
// L-scheme stability constant is updated by PKs.
// -----------------------------------------------------------------------------
template<class PK_Base>
std::vector<Key>
PK_MPC<PK_Base>::SetupLSchemeKey(Teuchos::ParameterList& plist)
{
  if (L_scheme_) {
    for (auto pk = sub_pks_.begin(); pk != sub_pks_.end(); ++pk) {
      auto tmp = (*pk)->SetupLSchemeKey(plist);
      L_scheme_keys_.insert(L_scheme_keys_.end(), tmp.begin(), tmp.end());
    }
  }
  return L_scheme_keys_;
}


// -----------------------------------------------------------------------------
// L-scheme initialization at the begging of time step
// -----------------------------------------------------------------------------
template<class PK_Base>
void
PK_MPC<PK_Base>::InitializeLSchemeStep()
{
  if (L_scheme_) {
    for (int i = 0; i < sub_pks_.size(); ++i) {
      if (!L_scheme_keys_[i].empty()) {
        Key key_prev = L_scheme_keys_[i] + "_prev";
        const auto& u0_c = *solution_->SubVector(i)->Data()->ViewComponent("cell");
        *S_->GetW<CompositeVector>(key_prev, "state").ViewComponent("cell") = u0_c;
      }
    }
  }
  auto& data = S_->GetW<LSchemeData>("l_scheme_data", "state");
  for (auto& item : data) item.second.reset();
}


// -----------------------------------------------------------------------------
// L-scheme stability constant is updated by PKs.
// -----------------------------------------------------------------------------
template<class PK_Base>
void
PK_MPC<PK_Base>::ComputeLSchemeStability()
{ 
  if (L_scheme_) {
    for (auto pk = sub_pks_.begin(); pk != sub_pks_.end(); ++pk) {
      (*pk)->ComputeLSchemeStability();
    }
  }
}

} // namespace Amanzi

#endif
