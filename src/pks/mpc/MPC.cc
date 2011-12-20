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
#include "MPC.hh"
#include "State.hh"


namespace Amanzi {

MPC::MPC(Teuchos::ParameterList &mpc_plist,
         Teuchos::RCP<State>& S, Teuchos::RCP<TreeVector>& soln) :
    mpc_plist_(mpc_plist), S_(S) {

  set_name(mpc_plist.get<std::string>("Name"));

  Teuchos::ParameterList pks_list = mpc_plist.sublist("PKs");
  for (Teuchos::ParameterList::ConstIterator i = pks_list.begin();
       i != pks_list.end(); ++i) {

    const std::string &name_i  = pks_list.name(i);
    const Teuchos::ParameterEntry  &entry_i = pks_list.entry(i);

    if (entry_i.isList()) {
      Teuchos::RCP<TreeVector> pk_soln = Teuchos::rcp(new TreeVector(name()));
      soln->PushBack(pk_soln);
      sub_pks_.push_back(pk_factory_.create_pk(pks_list.sublist(name_i),
                                               S, pk_soln));
    }
  }
};

void MPC::initialize(Teuchos::RCP<State>& S, Teuchos::RCP<TreeVector>& soln) {
  for (std::vector< Teuchos::RCP<PK> >::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {
    Tecuhos::RCP<TreeVector> pk_soln;
    int subvec_not_found = soln->SubVector((*pk)->name(), pk_soln);
    if (subvec_not_found) {
      Errors::Message message("MPC: vector structure does not match PK structure");
      Exceptions::amanzi_throw(message);
    }
    (*pk)->initialize(S, pk_soln);
  }
};

void MPC::state_to_solution(Teuchos::RCP<const State>& S,
                            Teuchos::RCP<TreeVector>& soln) {
  for (std::vector< Teuchos::RCP<PK> >::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {
    Tecuhos::RCP<TreeVector> pk_soln;
    int subvec_not_found = soln->SubVector((*pk)->name(), pk_soln);
    if (subvec_not_found) {
      Errors::Message message("MPC: vector structure does not match PK structure");
      Exceptions::amanzi_throw(message);
    }
    (*pk)->state_to_solution(S, pk_soln);
  }
};

void MPC::state_to_solution(Teuchos::RCP<const State>& S,
                            Teuchos::RCP<TreeVector>& soln,
                            Teuchos::RCP<TreeVector>& soln_dot) {
  for (std::vector< Teuchos::RCP<PK> >::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {
    Tecuhos::RCP<TreeVector> pk_soln;
    int subvec_not_found = soln->SubVector((*pk)->name(), pk_soln);
    if (subvec_not_found) {
      Errors::Message message("MPC: vector structure does not match PK structure");
      Exceptions::amanzi_throw(message);
    }

    Tecuhos::RCP<TreeVector> pk_soln_dot;
    int subvec_not_found = soln_dot->SubVector((*pk)->name(), pk_soln_dot);
    if (subvec_not_found) {
      Errors::Message message("MPC: vector structure does not match PK structure");
      Exceptions::amanzi_throw(message);
    }
    (*pk)->state_to_solution(S, pk_soln, pk_soln_dot);
  }
};

void MPC::solution_to_state(Teuchos::RCP<const TreeVector>& soln,
                            Teuchos::RCP<State>& S) {
  for (std::vector< Teuchos::RCP<PK> >::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {
    Tecuhos::RCP<const TreeVector> pk_soln;
    int subvec_not_found = soln->SubVector((*pk)->name(), pk_soln);
    if (subvec_not_found) {
      Errors::Message message("MPC: vector structure does not match PK structure");
      Exceptions::amanzi_throw(message);
    }
    (*pk)->solution_to_state(pk_soln, S);
  }
};

void MPC::solution_to_state(Teuchos::RCP<const TreeVector>& soln,
                            Teuchos::RCP<const TreeVector>& soln_dot,
                            Teuchos::RCP<State>& S) {
  for (std::vector< Teuchos::RCP<PK> >::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {
    Tecuhos::RCP<const TreeVector> pk_soln;
    int subvec_not_found = soln->SubVector((*pk)->name(), pk_soln);
    if (subvec_not_found) {
      Errors::Message message("MPC: vector structure does not match PK structure");
      Exceptions::amanzi_throw(message);
    }

    Tecuhos::RCP<const TreeVector> pk_soln_dot;
    int subvec_not_found = soln_dot->SubVector((*pk)->name(), pk_soln_dot);
    if (subvec_not_found) {
      Errors::Message message("MPC: vector structure does not match PK structure");
      Exceptions::amanzi_throw(message);
    }
    (*pk)->solution_to_state(pk_soln, pk_soln_dot, S);
  }
};

void state_to_solution(Teuchos::RCP<const State>& S,
                                 Teuchos::RCP<TreeVector>& soln,
                                 Teuchos::RCP<TreeVector>& soln_dot);
  virtual void solution_to_state(Teuchos::RCP<const TreeVector>& soln,
                                 Teuchos::RCP<State>& S);
  virtual void solution_to_state(Teuchos::RCP<const TreeVector>& soln,
                                 Teuchos::RCP<const TreeVector>& soln_dot,
                                 Teuchos::RCP<State>& S);

double MPC::get_dT() {
  double dt = 1.e99;
  for (std::vector< Teuchos::RCP<PK> >::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {
    dt = std::min(dt, (*pk)->get_dT());
  }
  return dt;
}

void MPC::compute_f(const double t, const Vector& u,
                    const Vector& udot, Vector& f) {
  const TreeVector *u_ptr = dynamic_cast<const TreeVector *> (&u);
  ASSERT(u_ptr);
  const TreeVector *udot_ptr = dynamic_cast<const TreeVector *> (&udot);
  ASSERT(udot_ptr);
  TreeVector *f_ptr = dynamic_cast<TreeVector *> (&f);
  ASSERT(f_ptr);

  // note this must happen outside of the pk loop -- it places the solution
  // for one PK in the state so that it may be accessed by other PKs
  solution_to_state(*u_ptr, *udot_ptr, S_);

  for (std::vector< Teuchos::RCP<PK> >::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {

    int subvec_not_found;
    Teuchos::RCP<const TreeVector> sub_u;
    subvec_not_found = u_ptr->SubVector((*pk)->name(), sub_u);
    if (subvec_not_found) {
      Errors::Message message("MPC: vector structure does not match PK structure");
      Exceptions::amanzi_throw(message);
    } 

    Teuchos::RCP<const TreeVector> sub_udot;
    subvec_not_found = udot_ptr->SubVector((*pk)->name(), sub_udot);
    if (subvec_not_found) {
      Errors::Message message("MPC: vector structure does not match PK structure");
      Exceptions::amanzi_throw(message);
    } 
    
    Teuchos::RCP<TreeVector> sub_f;
    subvec_not_found = f_ptr->SubVector((*pk)->name(), sub_f);
    if (subvec_not_found) {
      Errors::Message message("MPC: vector structure does not match PK structure");
      Exceptions::amanzi_throw(message);
    } 

    (*pk)->compute_f(t, *sub_u, *sub_udot, *sub_f);
  }
};

void MPC::commit_state(double dt, Teuchos::RCP<State> &S) {
  for (std::vector< Teuchos::RCP<PK> >::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {
    (*pk)->commit_state(dt, S);
  }
};
} // namespace
