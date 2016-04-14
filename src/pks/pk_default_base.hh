/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Default base with a few methods implemented in standard ways.
------------------------------------------------------------------------- */

#ifndef AMANZI_PK_DEFAULT_BASE_HH_
#define AMANZI_PK_DEFAULT_BASE_HH_

#include "Teuchos_ParameterList.hpp"

#include "VerboseObject.hh"
#include "PK.hh"


namespace Amanzi {

class PKDefaultBase : public PK {

 public:

  PKDefaultBase(const Teuchos::RCP<Teuchos::ParameterList>& plist,
                Teuchos::ParameterList& FElist,
                const Teuchos::RCP<TreeVector>& solution) :
      plist_(plist), solution_(solution) {}

  // Virtual destructor
  virtual ~PKDefaultBase() {}

  virtual void setup(const Teuchos::Ptr<State>& S);

  virtual void set_states(const Teuchos::RCP<const State>& S,
                          const Teuchos::RCP<State>& S_inter,
                          const Teuchos::RCP<State>& S_next);

  // -- ensure a solution is valid
  virtual bool valid_step() { return true; }

  virtual void solution_to_state(TreeVector& soln,
          const Teuchos::RCP<State>& S) = 0; // this is here to pass the buck virtually
  virtual void solution_to_state(const TreeVector& soln,
          const Teuchos::RCP<State>& S);


  virtual std::string name() { return name_; }

 protected:

  Teuchos::RCP<Teuchos::ParameterList> plist_;
  Teuchos::RCP<TreeVector> solution_;
  std::string name_;

  // states
  Teuchos::RCP<const State> S_;
  Teuchos::RCP<State> S_inter_;
  Teuchos::RCP<State> S_next_;

  // fancy OS
  Teuchos::RCP<VerboseObject> vo_;
};

} // namespace


#endif
