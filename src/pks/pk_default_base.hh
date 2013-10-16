/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Default base with a few methods implemented in standard ways.
------------------------------------------------------------------------- */

#ifndef AMANZI_PK_DEFAULT_BASE_HH_
#define AMANZI_PK_DEFAULT_BASE_HH_

#include "Teuchos_VerboseObject.hpp" // REMOVE ME
#include "Teuchos_ParameterList.hpp"

#include "VerboseObject.hh"
#include "PK.hh"


namespace Amanzi {

class PKDefaultBase : public PK {

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
  // REMOVE ME
  // bool includesVerbLevel(Teuchos::EVerbosityLevel my_level,
  //                        Teuchos::EVerbositylevel level, bool def) {
  //   return vo_->includesVerbLevel(my_level, level, def);
  // }

  // REMOVE ME!
  Teuchos::OSTab getOSTab(const int tabs=1) { return vo_->getOSTab(tabs); }

 protected:

  Teuchos::ParameterList plist_;
  Teuchos::RCP<TreeVector> solution_;
  std::string name_;

  // states
  Teuchos::RCP<const State> S_;
  Teuchos::RCP<State> S_inter_;
  Teuchos::RCP<State> S_next_;

  // fancy OS
  Teuchos::RCP<VerboseObject> vo_;

  // -- more fancy OS... these will go away eventually, but removal requires
  // -- much code refactoring.
  Teuchos::RCP<Teuchos::FancyOStream> out_;
  Teuchos::EVerbosityLevel verbosity_;

  // cruft for easier global debugging
  std::vector<AmanziMesh::Entity_ID> dc_;
  std::vector<Teuchos::RCP<VerboseObject> > dcvo_;
};

} // namespace


#endif
