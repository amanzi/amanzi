/* -*-  mode: c++; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

A DomainSet coupler, couples a bunch of domains of the same structure.

------------------------------------------------------------------------- */

#ifndef PKS_DOMAIN_SET_MPC_HH_
#define PKS_DOMAIN_SET_MPC_HH_

#include "Key.hh"
#include "PK.hh"
#include "mpc.hh"

namespace Amanzi {

class DomainSetMPC : public MPC<PK> {

 public:

  DomainSetMPC(Teuchos::ParameterList& pk_tree,
                           const Teuchos::RCP<Teuchos::ParameterList>& global_list,
                           const Teuchos::RCP<State>& S,
                           const Teuchos::RCP<TreeVector>& solution)
      : MPC<PK>(pk_tree, global_list, S, solution),
        PK(pk_tree, global_list, S, solution)
   {
    // grab the list of subpks
    auto subpks = this->plist_->template get<Teuchos::Array<std::string> >("PKs order");
    std::string dsname = subpks.back();
    subpks.pop_back();

    KeyTriple triple;
    bool is_ds = Keys::splitDomainSet(dsname, triple);
    if (!is_ds || !S->HasDomainSet(std::get<0>(triple))) {
      Errors::Message msg;
      msg << "DomainSetMPC: \"" << dsname << "\" should be a domain-set of the form DOMAIN_*-PK_NAME";
      Exceptions::amanzi_throw(msg);
    }

    // add for the various sub-pks based on IDs
    auto ds = S->GetDomainSet(std::get<0>(triple));
    for (auto& name_id : *ds) {
      subpks.push_back(Keys::getKey(name_id.first, std::get<2>(triple)));
    }
    this->plist_->template set("PKs order", subpks);

    // construct the sub-PKs on COMM_SELF
    MPC<PK>::init_(S, getCommSelf());
  }

  // Virtual destructor
  virtual ~DomainSetMPC() = default;
  
  // PK methods
  // -- dt is the minimum of the sub pks
  virtual double get_dt();

  // -- advance each sub pk dt.
  virtual bool AdvanceStep(double t_old, double t_new, bool reinit);

  virtual void set_dt(double dt);


 protected:
  std::string pks_set_;

 private:
  // factory registration
  static RegisteredPKFactory<DomainSetMPC> reg_;

};

} // namspace
        
#endif
