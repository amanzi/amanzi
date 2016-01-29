/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------


Author: Ethan Coon

Default base with a few methods implemented in standard ways.
------------------------------------------------------------------------- */

#ifndef AMANZI_PK_DEFAULT_BASE_HH_
#define AMANZI_PK_DEFAULT_BASE_HH_

#include "Teuchos_ParameterList.hpp"

#include "VerboseObject.hh"
#include "PK.hh"

namespace Amanzi {

class TreeVector;

class PKDefaultBase : public PK {

public:

  PKDefaultBase() {};
  PKDefaultBase(Teuchos::ParameterList& pk_tree,
                 const Teuchos::RCP<Teuchos::ParameterList>& glist,
                 const Teuchos::RCP<State>& S,
                 const Teuchos::RCP<TreeVector>& soln) {};

  // Virtual destructor
  virtual ~PKDefaultBase() {}

  virtual void Setup() {}

  virtual void CommitStep(double t_old, double t_new) {}

  virtual void CalculateDiagnostics(){}

  // virtual void set_states(const Teuchos::RCP<const State>& S,
  //                         const Teuchos::RCP<State>& S_inter,
  //                         const Teuchos::RCP<State>& S_next);


  // virtual void solution_to_state(TreeVector& soln,
  //         const Teuchos::RCP<State>& S) = 0; // this is here to pass the buck virtually
  // virtual void solution_to_state(const TreeVector& soln,
  //         const Teuchos::RCP<State>& S);


  //virtual std::string name() { return name_; }
  virtual std::string name() {}

protected:

 //  Teuchos::RCP<Teuchos::ParameterList> plist_;
 //  Teuchos::RCP<TreeVector> solution_;
 //  std::string name_;

 //  // states
 //  Teuchos::RCP<const State> S_;
 //  Teuchos::RCP<State> S_inter_;
 //  Teuchos::RCP<State> S_next_;

 //  // fancy OS
  VerboseObject* vo_;
};

} // namespace


#endif
