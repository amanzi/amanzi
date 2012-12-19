/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Default base with a few methods implemented in standard ways.
------------------------------------------------------------------------- */

#ifndef AMANZI_PK_DEFAULT_BASE_HH_
#define AMANZI_PK_DEFAULT_BASE_HH_

#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"
#include "Teuchos_ParameterList.hpp"
#include "PK.hh"


namespace Amanzi {

class PKDefaultBase : public PK, public Teuchos::VerboseObject<PKDefaultBase> {

 public:

  PKDefaultBase(Teuchos::ParameterList& plist,
                const Teuchos::RCP<TreeVector>& solution) :
      plist_(plist), solution_(solution) {}

  // Virtual destructor
  virtual ~PKDefaultBase() {}

  virtual void setup(const Teuchos::Ptr<State>& S);

  virtual void set_states(const Teuchos::RCP<const State>& S,
                          const Teuchos::RCP<State>& S_inter,
                          const Teuchos::RCP<State>& S_next);


  virtual void solution_to_state(const Teuchos::RCP<TreeVector>& soln,
          const Teuchos::RCP<State>& S) = 0; // this is here to pass the buck virtually
  virtual void solution_to_state(const Teuchos::RCP<const TreeVector>& soln,
          const Teuchos::RCP<State>& S);


  virtual std::string name() { return name_; }

 protected:

  Teuchos::ParameterList plist_;
  Teuchos::RCP<TreeVector> solution_;
  std::string name_;

  // states
  Teuchos::RCP<const State> S_;
  Teuchos::RCP<State> S_inter_;
  Teuchos::RCP<State> S_next_;

  // fancy OS
  Teuchos::RCP<Teuchos::FancyOStream> out_;
  Teuchos::EVerbosityLevel verbosity_;
};

} // namespace


#endif
