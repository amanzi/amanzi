/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Default base with a few methods implemented in standard ways.
------------------------------------------------------------------------- */
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "pk_default_base.hh"

namespace Amanzi {


void PKDefaultBase::setup(const Teuchos::Ptr<State>& S) {
  // THIS MAY BE CALLED MORE THAN ONCE!
  name_ = plist_.get<std::string>("PK name");

  // set up the VerboseObject
  vo_ = Teuchos::rcp(new VerboseObject(name_, plist_));

  // Mirror the vo_ in the pk, as removal will require much code refactoring
  out_ = vo_->os();
  verbosity_ = vo_->getVerbLevel();


  // cells to debug
  if (dc_.size() == 0) {
    if (plist_.isParameter("debug cells")) {
      Teuchos::Array<int> dc = plist_.get<Teuchos::Array<int> >("debug cells");
      for (Teuchos::Array<int>::const_iterator c=dc.begin();
           c!=dc.end(); ++c) dc_.push_back(*c);

      // Enable a vo for each cell, allows parallel printing of debug cells.
      if (plist_.isParameter("debug cell ranks")) {
        Teuchos::Array<int> dc_ranks = plist_.get<Teuchos::Array<int> >("debug cell ranks");
        if (dc.size() != dc_ranks.size()) {
          Errors::Message message("Debug cell and debug cell ranks must be equal length.");
          Exceptions::amanzi_throw(message);
        }
        for (Teuchos::Array<int>::const_iterator dcr=dc_ranks.begin();
             dcr!=dc_ranks.end(); ++dcr) {
          // make a verbose object for each case
          Teuchos::ParameterList vo_plist;
          vo_plist.sublist("VerboseObject");
          vo_plist.sublist("VerboseObject")
              = plist_.sublist("VerboseObject");
          vo_plist.sublist("VerboseObject").set("Root ID", *dcr);

          Teuchos::writeParameterListToXmlOStream(vo_plist, std::cout);
          dcvo_.push_back(Teuchos::rcp(new VerboseObject(name_,vo_plist)));

        }
      } else {
        // Simply use the pk's vo
        dcvo_.resize(dc_.size(), vo_);
      }
    }
  }

}


void PKDefaultBase::set_states(const Teuchos::RCP<const State>& S,
        const Teuchos::RCP<State>& S_inter,
        const Teuchos::RCP<State>& S_next) {
  S_ = S;
  S_inter_ = S_inter;
  S_next_ = S_next;
}


void PKDefaultBase::solution_to_state(const Teuchos::RCP<const TreeVector>& soln,
        const Teuchos::RCP<State>& S) {
  Teuchos::RCP<TreeVector> nc_soln = Teuchos::rcp_const_cast<TreeVector>(soln);
  solution_to_state(nc_soln, S);
}

} // namespace
