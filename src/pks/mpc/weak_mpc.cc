/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Implementation for the derived WeakMPC class.  Provides only the advance()
method missing from MPC.hh.  In weak coupling, we simply loop over the
sub-PKs, calling their advance() methods and returning failure if any fail.

See additional documentation in the base class src/pks/mpc/MPC.hh
------------------------------------------------------------------------- */

#include <string>

#include "errors.hh"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"

#include "weak_mpc.hh"

namespace Amanzi {

RegisteredPKFactory<WeakMPC> WeakMPC::reg_("weak MPC");

WeakMPC::WeakMPC(Teuchos::ParameterList &mpc_plist, const Teuchos::RCP<State>& S,
         const Teuchos::RCP<TreeVector>& solution) :
    MPC(mpc_plist, S, solution) {

  solution_ = solution;

  // loop over sub-PKs in the PK sublist, constructing the hierarchy recursively
  Teuchos::ParameterList pks_list = mpc_plist.sublist("PKs");
  PKFactory pk_factory;

  if (mpc_plist.isParameter("PKs order")) {
    // ordered
    Teuchos::Array<std::string> pk_order = mpc_plist.get< Teuchos::Array<std::string> >("PKs order");
    unsigned int npks = pk_order.size();

    for (unsigned int i=0; i!=npks; ++i) {
      std::string name_i = pk_order[i];
      Teuchos::RCP<TreeVector> pk_soln = Teuchos::rcp(new TreeVector(name_i));
      solution_->PushBack(pk_soln);
      Teuchos::RCP<PK> pk = pk_factory.CreatePK(pks_list.sublist(name_i), S, pk_soln);
      pk->set_name(name_i);
      sub_pks_.push_back(pk);
    }

  } else {
    // no order, just loop over all sublists
    for (Teuchos::ParameterList::ConstIterator i = pks_list.begin();
         i != pks_list.end(); ++i) {

      const std::string &name_i  = pks_list.name(i);
      const Teuchos::ParameterEntry  &entry_i = pks_list.entry(i);
      if (entry_i.isList()) {
        Teuchos::RCP<TreeVector> pk_soln = Teuchos::rcp(new TreeVector(name_i));
        solution_->PushBack(pk_soln);

        Teuchos::RCP<PK> pk = pk_factory.CreatePK(pks_list.sublist(name_i), S, pk_soln);
        pk->set_name(name_i);
        sub_pks_.push_back(pk);
      }
    }
  }
};


// Advance each sub-PK individually.
bool WeakMPC::advance(double dt) {
  bool fail = false;
  for (std::vector< Teuchos::RCP<PK> >::iterator pk = sub_pks_.begin();
       pk != sub_pks_.end(); ++pk) {
    fail = (*pk)->advance(dt);
    if (fail) {
      return fail;
    }
  }
  return fail;
};
} // namespace
