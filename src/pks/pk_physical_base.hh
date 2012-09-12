/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Default base with default implementations of methods for a physical PK.
------------------------------------------------------------------------- */

#ifndef AMANZI_PK_PHYSICAL_BASE_HH_
#define AMANZI_PK_PHYSICAL_BASE_HH_

#include "primary_variable_field_evaluator.hh"
#include "pk_default_base.hh"

namespace Amanzi {

class PKPhysicalBase : public virtual PKDefaultBase {

 public:
  PKPhysicalBase(Teuchos::ParameterList& plist,
                 const Teuchos::RCP<TreeVector>& solution) :
      PKDefaultBase(plist,solution) {}

  // Default implementations of PK methods.
  // -- transfer operators -- pointer copies only
  virtual void state_to_solution(const Teuchos::RCP<State>& S,
                                 const Teuchos::RCP<TreeVector>& soln);
  virtual void solution_to_state(const Teuchos::RCP<TreeVector>& soln,
                                 const Teuchos::RCP<State>& S);


  // -- setup
  virtual void setup(const Teuchos::Ptr<State>& S);

  // -- initialize
  virtual void initialize(const Teuchos::Ptr<State>& S);

 protected: // methods
  std::string Key_(std::string suffix) { return domain_prefix_+suffix; }

 protected: // data
  // name of domain, associated mesh
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  std::string domain_;
  std::string domain_prefix_;

  // solution and evaluator
  std::string key_;
  Teuchos::RCP<PrimaryVariableFieldEvaluator> solution_evaluator_;

};

} // namespace

#endif
