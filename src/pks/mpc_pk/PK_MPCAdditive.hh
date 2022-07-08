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

#ifndef AMANZI_MPCADD_HH_
#define AMANZI_MPCADD_HH_

#include <string>
#include <vector>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Epetra_MpiComm.h"

#include "State.hh"
#include "TreeVector.hh"

#include "PK_MPC.hh"
#include "PK_Factory.hh"

namespace Amanzi {

template <class PK_Base>
class PK_MPCAdditive :  public PK_MPC<PK> {
 public:
  PK_MPCAdditive(Teuchos::ParameterList& pk_tree,
         const Teuchos::RCP<Teuchos::ParameterList>& global_list,
         const Teuchos::RCP<State>& S,
         const Teuchos::RCP<TreeVector>& soln);

  ~PK_MPCAdditive() {};

  // PK methods
  // -- sets up sub-PKs
  virtual void Setup();

  // -- calls all sub-PK initialize() methods
  virtual void Initialize();

  // -- loops over sub-PKs
  virtual void CommitStep(double t_old, double t_new, const Tag& tag);
  virtual void CalculateDiagnostics(const Tag& tag);
  virtual double get_dt();
  virtual void set_dt(double dt);
  virtual bool AdvanceStep(double t_old, double t_new, bool reinit = false);
};


// -----------------------------------------------------------------------------
// Setup of PK hierarchy from PList
// -----------------------------------------------------------------------------
template <class PK_Base>
PK_MPCAdditive<PK_Base>::PK_MPCAdditive(
    Teuchos::ParameterList& pk_tree,
    const Teuchos::RCP<Teuchos::ParameterList>& global_list,
    const Teuchos::RCP<State>& S,
    const Teuchos::RCP<TreeVector>& soln)
{
  pk_tree_ = pk_tree;
  global_list_ = global_list;
  S_ = S;
  solution_ = soln;

  // name the PK
  name_ = pk_tree.name();
  auto found = name_.rfind("->");
  if (found != std::string::npos) name_.erase(0, found + 2);

  // const char* result = name_.data();
  // while ((result = std::strstr(result, "->")) != NULL) {
  //   result += 2;
  //   name_ = result;   
  // }

  // get my ParameterList
  my_list_ = Teuchos::sublist(Teuchos::sublist(global_list_, "PKs"), name_);

  Teuchos::RCP<Teuchos::ParameterList> plist;
  if (global_list_->isSublist("PKs")) {
    plist = Teuchos::sublist(global_list, "PKs");
  }

  std::vector<std::string> pk_name = my_list_->get<Teuchos::Array<std::string> >("PKs order").toVector();

  // loop over sub-PKs in the PK sublist, constructing the hierarchy recursively
  PKFactory pk_factory;

  // Create only one subvector for all PKs
  Teuchos::RCP<TreeVector> pk_soln = Teuchos::rcp(new TreeVector());
  solution_->PushBack(pk_soln);

  sub_pks_.clear();
  for (int i = 0; i < pk_name.size(); i++) {
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
void PK_MPCAdditive<PK_Base>::Setup() {
  for (typename SubPKList::iterator pk = sub_pks_.begin(); pk != sub_pks_.end(); ++pk) {
    (*pk)->Setup();
  }
}


// -----------------------------------------------------------------------------
// Loop over sub-PKs, calling their initialization methods
// -----------------------------------------------------------------------------
template <class PK_Base>
void PK_MPCAdditive<PK_Base>::Initialize() {
  for (typename SubPKList::iterator pk = sub_pks_.begin(); pk != sub_pks_.end(); ++pk) {
    (*pk)->Initialize();
  }
}


// -----------------------------------------------------------------------------
// loop over sub-PKs, calling their commit state method
// -----------------------------------------------------------------------------
template <class PK_Base>
void PK_MPCAdditive<PK_Base>::CommitStep(double t_old, double t_new, const Tag& tag) {
  for (typename SubPKList::iterator pk = sub_pks_.begin(); pk != sub_pks_.end(); ++pk) {
    (*pk)->CommitStep(t_old, t_new, tag);
  }
}


// -----------------------------------------------------------------------------
// loop over sub-PKs, calling their CalculateDiagnostics method
// -----------------------------------------------------------------------------
template <class PK_Base>
void PK_MPCAdditive<PK_Base>::CalculateDiagnostics(const Tag& tag) {
  for (typename SubPKList::iterator pk = sub_pks_.begin(); pk != sub_pks_.end(); ++pk) {
    (*pk)->CalculateDiagnostics(tag);
  }
}


// -----------------------------------------------------------------------------
// Calculate the min of sub PKs timestep sizes.
// -----------------------------------------------------------------------------
template <class PK_Base>
double PK_MPCAdditive<PK_Base>::get_dt() {
  double dt = 1.0e99;
  for (PK_MPCAdditive<PK>::SubPKList::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {
    dt = std::min<double>(dt, (*pk)->get_dt());
  }
  return dt;
}

// -----------------------------------------------------------------------------
// Set sub PKs timestep sizes.
// -----------------------------------------------------------------------------
template <class PK_Base>
void PK_MPCAdditive<PK_Base>::set_dt(double dt_) {
  // double dt = 1.0e99;
  for (PK_MPCAdditive<PK>::SubPKList::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {
    (*pk)->set_dt(dt_);
  }
}


// -----------------------------------------------------------------------------
// Advance each sub-PK individually, returning a failure as soon as possible.
// -----------------------------------------------------------------------------
template <class PK_Base>
bool PK_MPCAdditive<PK_Base>::AdvanceStep(double t_old, double t_new, bool reinit) {

  bool fail = false;
  for (PK_MPCAdditive<PK>::SubPKList::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {
    fail = (*pk)->AdvanceStep(t_old, t_new, reinit);
    if (fail) {
      return fail;
    }
  }
  return fail;

}


}  // namespace Amanzi

#endif
