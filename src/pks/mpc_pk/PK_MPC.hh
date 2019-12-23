/*
  This is the mpc_pk component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
           Daniil Svyatskiy

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

#include "State.hh"
#include "TreeVector.hh"

#include "PK.hh"
#include "PK_Factory.hh"

namespace Amanzi {

template <class PK_Base>
class PK_MPC : virtual public PK {
 public:
  PK_MPC(Teuchos::ParameterList& pk_tree,
         const Teuchos::RCP<Teuchos::ParameterList>& global_list,
         const Teuchos::RCP<State>& S,
         const Teuchos::RCP<TreeVector>& soln);

  PK_MPC() {};

  ~PK_MPC() {};

  // PK methods
  // -- sets up sub-PKs
  virtual void Setup(const Teuchos::Ptr<State>& S);

  // -- calls all sub-PK initialize() methods
  virtual void Initialize(const Teuchos::Ptr<State>& S);

  // -- loops over sub-PKs
  virtual void CommitStep(double t_old, double t_new, const Teuchos::RCP<State>& S);
  virtual void CalculateDiagnostics(const Teuchos::RCP<State>& S);

  virtual void set_states(const Teuchos::RCP<const State>& S,
                          const Teuchos::RCP<State>& S_inter,
                          const Teuchos::RCP<State>& S_next);

  virtual void Solution_to_State(TreeVector& soln,
                                 const Teuchos::RCP<State>& S) {};
  virtual void Solution_to_State(const TreeVector& soln,
                                 const Teuchos::RCP<State>& S) {};

  virtual void State_to_Solution(const Teuchos::RCP<State>& S,
                                 TreeVector& soln) {} ;

  int num_sub_pks() {return sub_pks_.size();}
  // -- identifier accessor
  std::string name() const { return name_; }
  virtual Teuchos::RCP<PK_Base> get_subpk(int i);
    
  
 protected:
  // identifier
  std::string name_;

  // list of the PKs coupled by this MPC
  typedef std::vector<Teuchos::RCP<PK_Base> > SubPKList;
  SubPKList sub_pks_;

  // single solution vector for the coupled problem
  Teuchos::RCP<TreeVector> solution_;

  // plists
  Teuchos::RCP<Teuchos::ParameterList> global_list_;
  Teuchos::RCP<Teuchos::ParameterList> my_list_;
  Teuchos::ParameterList pk_tree_;

  // states
  Teuchos::RCP<State> S_;
};


// -----------------------------------------------------------------------------
// Setup of PK hierarchy from PList
// -----------------------------------------------------------------------------
template <class PK_Base>
PK_MPC<PK_Base>::PK_MPC(Teuchos::ParameterList& pk_tree,
                        const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                        const Teuchos::RCP<State>& S,
                        const Teuchos::RCP<TreeVector>& soln) :
  pk_tree_(pk_tree),
  global_list_(global_list),
  S_(S),
  solution_(soln)
{
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

  std::vector<std::string> pk_name = my_list_->get<Teuchos::Array<std::string> >("PKs order").toVector();

  // loop over sub-PKs in the PK sublist, constructing the hierarchy recursively
  PKFactory pk_factory;

  for (int i = 0; i < pk_name.size(); ++i) {
    //const std::string& sub_name = sub->first;
    const std::string& sub_name = pk_name[i];
    if (!plist->isSublist(sub_name)) {
      Errors::Message msg;
      msg << "PK tree has no sublist \"" << sub_name << "\"";
      Exceptions::amanzi_throw(msg);
    }
  }

  for (int i = 0; i < pk_name.size(); i++) {
    // Collect arguments to the constructor
    Teuchos::RCP<TreeVector> pk_soln = Teuchos::rcp(new TreeVector());
    solution_->PushBack(pk_soln);

    // create the PK
    Teuchos::RCP<PK> pk_notype = pk_factory.CreatePK(pk_name[i], pk_tree, global_list, S, pk_soln);
    Teuchos::RCP<PK_Base> pk = Teuchos::rcp_dynamic_cast<PK_Base>(pk_notype);
    sub_pks_.push_back(pk);
  }
}


// -----------------------------------------------------------------------------
// Setup of PK hierarchy from PList
// -----------------------------------------------------------------------------
template <class PK_Base>
void PK_MPC<PK_Base>::Setup(const Teuchos::Ptr<State>& S) {
  for (typename SubPKList::iterator pk = sub_pks_.begin(); pk != sub_pks_.end(); ++pk) {
    (*pk)->Setup(S);
  }
}


// -----------------------------------------------------------------------------
// Loop over sub-PKs, calling their initialization methods
// -----------------------------------------------------------------------------
template <class PK_Base>
void PK_MPC<PK_Base>::Initialize(const Teuchos::Ptr<State>& S) {
  for (typename SubPKList::iterator pk = sub_pks_.begin(); pk != sub_pks_.end(); ++pk) {
    (*pk)->Initialize(S);
  }
}


// -----------------------------------------------------------------------------
// loop over sub-PKs, calling their commit state method
// -----------------------------------------------------------------------------
template <class PK_Base>
void PK_MPC<PK_Base>::CommitStep(double t_old, double t_new, const Teuchos::RCP<State>& S) {
  for (typename SubPKList::iterator pk = sub_pks_.begin(); pk != sub_pks_.end(); ++pk) {
    (*pk)->CommitStep(t_old, t_new, S);
  }
}


// -----------------------------------------------------------------------------
// loop over sub-PKs, calling their CalculateDiagnostics method
// -----------------------------------------------------------------------------
template <class PK_Base>
void PK_MPC<PK_Base>::CalculateDiagnostics(const Teuchos::RCP<State>& S) {
  for (typename SubPKList::iterator pk = sub_pks_.begin(); pk != sub_pks_.end(); ++pk) {
    (*pk)->CalculateDiagnostics(S);
  }
}

template <class PK_Base>
void PK_MPC<PK_Base> :: set_states(const Teuchos::RCP<const State>& S,
                                   const Teuchos::RCP<State>& S_inter,
                                   const Teuchos::RCP<State>& S_next) {

  for (typename SubPKList::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {
    (*pk)->set_states(S, S_inter, S_next);
  }

};


template <class PK_Base>
Teuchos::RCP<PK_Base> PK_MPC<PK_Base> :: get_subpk(int i){

  if (i >= sub_pks_.size()) {
    return Teuchos::null;
  }else{
    return sub_pks_.at(i);
  }

}  

}  // namespace Amanzi

#endif
