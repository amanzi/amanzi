/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------


Author: Ethan Coon

Default base with a few methods implemented in standard ways.
------------------------------------------------------------------------- */

#ifndef AMANZI_PK_PHYSICAL_BASE_HH_
#define AMANZI_PK_PHYSICAL_BASE_HH_

#include "Teuchos_ParameterList.hpp"

#include "VerboseObject.hh"
#include "primary_variable_field_evaluator.hh"
#include "PK_default_base.hh"

namespace Amanzi {

class PKPhysicalBase : public virtual PKDefaultBase {

 public:

  PKPhysicalBase(){};

  PKPhysicalBase(Teuchos::ParameterList& pk_tree,
                 const Teuchos::RCP<Teuchos::ParameterList>& glist,
                 const Teuchos::RCP<State>& S,
                 const Teuchos::RCP<TreeVector>& soln){};


  // Virtual destructor
  virtual ~PKPhysicalBase() {}

  // Default implementations of PK methods.
  // -- transfer operators -- pointer copies only
  // virtual void state_to_solution(const Teuchos::RCP<State>& S,
  //                                TreeVector& soln);
  // virtual void solution_to_state(TreeVector& soln,
  //                                const Teuchos::RCP<State>& S);

  // new virtual set_states() to also get the primary field evaulator.
  // virtual void set_states(const Teuchos::RCP<const State>& S,
  //         const Teuchos::RCP<State>& S_inter,
  //         const Teuchos::RCP<State>& S_next);


  // -- setup
  virtual void Setup(){};

  // -- initialize
  virtual void Initialize(){};

  // Accessor for debugger, for use by coupling MPCs
  //  Teuchos::RCP<Debugger> debugger() { return db_; }

 protected: // methods
  //  void DeriveFaceValuesFromCellValues_(const Teuchos::Ptr<CompositeVector>& cv);

 protected: // data
  // name of domain, associated mesh
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  std::string domain_;

  // solution and evaluator
  std::string key_;
  //  Teuchos::RCP<PrimaryVariableFieldEvaluator> solution_evaluator_;

  // debugger for dumping vectors
  //Teuchos::RCP<Debugger> db_;

  // ENORM struct
  // typedef struct ENorm_t {
  //   double value;
  //   int gid;
  // } ENorm_t;

};

} // namespace

#endif
