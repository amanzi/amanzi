/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

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
#include "WeakMPC.hh"
#include "State.hh"


namespace Amanzi {

WeakMPC::WeakMPC(Teuchos::ParameterList &mpc_plist,
                 Teuchos::RCP<State> &S) :
    mpc_plist_(mpc_plist), S_(S) {

  set_name(mpc_plist.get<std::string>("Name"));

  Teuchos::ParameterList pks_list = mpc_plist.sublist("PKs");
  for (Teuchos::ParameterList::ConstIterator i = pks_list.begin();
       i != pks_list.end(); ++i) {

    const std::string &name_i  = pks_list.name(i);
    const Teuchos::ParameterEntry  &entry_i = pks_list.entry(i);

    if (entry_i.isList()) {
      sub_pks_.push_back(pk_factory_.create_pk(pks_list.sublist(name_i),S));
    }
  }
};

void WeakMPC::initialize() {
  for (std::vector< Teuchos::RCP<PK> >::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {
    (*pk)->initialize();
  }
};

double WeakMPC::get_dT() {
  double dt = 1.e99;
  for (std::vector< Teuchos::RCP<PK> >::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {
    dt = std::min(dt, (*pk)->get_dT());
  }
  return dt;
}

bool WeakMPC::advance_transient(double dt, Teuchos::RCP<State> &S0,
             Teuchos::RCP<State> &S1, Teuchos::RCP<TreeVector> &solution) {
  // check that the solution vector matches the sub_pks

  bool fail = 0;
  for (std::vector< Teuchos::RCP<PK> >::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {
    Techos::RCP<TreeVector> subvec;
    int pk_not_found = solution->SubVector((*pk)->name(), subvec);
    if (pk_not_found) {
      Errors::Message message("WeakMPC: vector structure does not match PK structure");
      Exceptions::amanzi_throw(message);
    } else {
      fail = (*pk)->advance_transient(dt, S0, S1, subvec);
      if (fail) {
        return fail;
      }
    }
  }
  return fail;
};

void WeakMPC::commit_state(double dt, Teuchos::RCP<State> &S) {
  for (std::vector< Teuchos::RCP<PK> >::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {
    (*pk)->commit_state(dt, S);
  }
};

void WeakMPC::solution_to_state(const TreeVector& u,
                                const TreeVector& udot) {
  for (std::vector< Teuchos::RCP<PK> >::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {
    int pk_not_found;

    Techos::RCP<TreeVector> sub_u;
    pk_not_found = u->SubVector((*pk)->name(), sub_u);
    if (pk_not_found) {
      Errors::Message message("WeakMPC: vector structure does not match PK structure");
      Exceptions::amanzi_throw(message);
    } 

    Techos::RCP<TreeVector> sub_udot;
    pk_not_found = udot->SubVector((*pk)->name(), sub_udot);
    if (pk_not_found) {
      Errors::Message message("WeakMPC: vector structure does not match PK structure");
      Exceptions::amanzi_throw(message);
    } 

    (*pk)->solution_to_state(*sub_u, *sub_udot);
  }
};

void WeakMPC::compute_f(const double t, const Vector& u,
                        const Vector& udot, Vector& f) {
  const TreeVector *u_ptr = dynamic_cast<const TreeVector *> (&u);
  ASSERT(u_ptr);
  const TreeVector *udot_ptr = dynamic_cast<const TreeVector *> (&udot);
  ASSERT(udot_ptr);
  TreeVector *f_ptr = dynamic_cast<TreeVector *> (&f);
  ASSERT(f_ptr);

  solution_to_state(*u_ptr, *udot_ptr);
  for (std::vector< Teuchos::RCP<PK> >::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {

    Techos::RCP<TreeVector> sub_u;
    pk_not_found = u_ptr->SubVector((*pk)->name(), *sub_u);
    if (pk_not_found) {
      Errors::Message message("WeakMPC: vector structure does not match PK structure");
      Exceptions::amanzi_throw(message);
    } 

    Techos::RCP<TreeVector> sub_udot;
    pk_not_found = udot_ptr->SubVector((*pk)->name(), *sub_udot);
    if (pk_not_found) {
      Errors::Message message("WeakMPC: vector structure does not match PK structure");
      Exceptions::amanzi_throw(message);
    } 
    
    Techos::RCP<TreeVector> sub_f;
    pk_not_found = f_ptr->SubVector((*pk)->name(), *sub_f);
    if (pk_not_found) {
      Errors::Message message("WeakMPC: vector structure does not match PK structure");
      Exceptions::amanzi_throw(message);
    } 

    (*pk)->compute_f(t, *sub_u, *sub_udot, *sub_f);

} // namespace
