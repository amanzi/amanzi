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
#include "State.hh"

#include "WeakMPC.hh"

namespace Amanzi {

WeakMPC::WeakMPC(Teuchos::ParameterList& mpc_plist,
                 Teuchos::RCP<State>& S) :
    mpc_plist_(mpc_plist), S_(S) {
  MPC::MPC(Teuchos::ParameterList& mpc_plist, Teuchos::RCP<State>& S) {
};

bool WeakMPC::advance(double dt, Teuchos::RCP<State>& S0,
             Teuchos::RCP<State>& S1, Teuchos::RCP<TreeVector>& solution) {
  bool fail = false;
  for (std::vector< Teuchos::RCP<PK> >::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {
    Teuchos::RCP<TreeVector> subvec;
    int subvec_not_found = solution->SubVector((*pk)->name(), subvec);
    if (subvec_not_found) {
      Errors::Message message("WeakMPC: vector structure does not match PK structure");
      Exceptions::amanzi_throw(message);
    } else {
      fail = (*pk)->advance(dt, S0, S1, subvec);
      if (fail) {
        return fail;
      }
    }
  }
  return fail;
};
} // namespace
