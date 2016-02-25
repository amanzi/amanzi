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

class PK_Default : public PK {

public:

  PK_Default() {};

  PK_Default(Teuchos::ParameterList& pk_tree,
             const Teuchos::RCP<Teuchos::ParameterList>& glist,
             const Teuchos::RCP<State>& S,
             const Teuchos::RCP<TreeVector>& soln):
    solution_(soln){};

  PK_Default(const Teuchos::RCP<Teuchos::ParameterList>& plist,
             Teuchos::ParameterList& FElist,
             const Teuchos::RCP<TreeVector>& solution) :
      plist_(plist), solution_(solution) {}

  // Virtual destructor
  virtual ~PK_Default() {}

  virtual void Setup(const Teuchos::Ptr<State>& S);

  virtual void set_states(const Teuchos::RCP<const State>& S,
                           const Teuchos::RCP<State>& S_inter,
                           const Teuchos::RCP<State>& S_next);


  // virtual void Solution_to_State(TreeVector& soln,
  //         const Teuchos::RCP<State>& S) = 0; 
  virtual void Solution_to_State(const TreeVector& soln,
          const Teuchos::RCP<State>& S);


  virtual std::string name() { return name_; }
  //virtual std::string name() {}

protected:

  Teuchos::RCP<Teuchos::ParameterList> plist_;
  Teuchos::RCP<TreeVector> solution_;
  std::string name_;

 //  // states
  Teuchos::RCP<const State> S_;
  Teuchos::RCP<State> S_inter_;
  Teuchos::RCP<State> S_next_;

 //  // fancy OS
  Teuchos::RCP<VerboseObject> vo_;
  //VerboseObject* vo_;
};

} // namespace


#endif
