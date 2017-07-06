/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------


Author: Ethan Coon

Default base with a few methods implemented in standard ways.
------------------------------------------------------------------------- */

#ifndef AMANZI_PK_PHYSICAL_HH_
#define AMANZI_PK_PHYSICAL_HH_

#include "Teuchos_ParameterList.hpp"
#include "VerboseObject.hh"
#include "primary_variable_field_evaluator.hh"
#include "PK.hh"
#include "Debugger.hh"

namespace Amanzi {

class PK_Physical : virtual public PK {
 public:
  PK_Physical() {}

  PK_Physical(Teuchos::ParameterList& pk_tree,
              const Teuchos::RCP<Teuchos::ParameterList>& glist,
              const Teuchos::RCP<State>& S,
              const Teuchos::RCP<TreeVector>& soln):
      PK(pk_tree, glist, S, soln) {};


  // Virtual destructor
  virtual ~PK_Physical() {}

  // Default implementations of PK methods.
  // -- transfer operators -- pointer copies only
  virtual void State_to_Solution(const Teuchos::RCP<State>& S,
                                  TreeVector& soln);
  virtual void Solution_to_State(TreeVector& soln,
                                  const Teuchos::RCP<State>& S);
  virtual void Solution_to_State(const TreeVector& soln,
                                  const Teuchos::RCP<State>& S);

  // new virtual set_states() to also get the primary field evaulator.
  virtual void set_states(const Teuchos::RCP<const State>& S,
                          const Teuchos::RCP<State>& S_inter,
                          const Teuchos::RCP<State>& S_next);

  virtual void ChangedSolutionPK() {
    solution_evaluator_->SetFieldAsChanged(S_next_.ptr());
  }
  
  // Accessor for debugger, for use by coupling MPCs
  Teuchos::RCP<Debugger> debugger() { return db_; }

 protected:
  // name of domain, associated mesh
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  std::string domain_;

  // solution and evaluator
  std::string key_;
  Teuchos::RCP<PrimaryVariableFieldEvaluator> solution_evaluator_;

  // debugger for dumping vectors
  Teuchos::RCP<Debugger> db_;
};

} // namespace Amanzi

#endif
