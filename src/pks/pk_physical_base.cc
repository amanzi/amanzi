/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
ATS

License: see $ATS_DIR/COPYRIGHT
Author: Ethan Coon

Default base with default implementations of methods for a physical PK.
------------------------------------------------------------------------- */

#include "pk_physical_base.hh"

namespace Amanzi {


// -----------------------------------------------------------------------------
// Construction of data.
// -----------------------------------------------------------------------------
void PKPhysicalBase::setup(const Teuchos::Ptr<State>& S) {
  PKDefaultBase::setup(S);

  // domain (defaults to nothing)
  domain_ = plist_.get<std::string>("domain name", std::string("domain"));
  if (domain_ != std::string("domain")) {
    domain_prefix_ = domain_ + std::string("_");
  }

  // get the mesh
  mesh_ = S->GetMesh(domain_);

  // solution key
  key_ = plist_.get<std::string>("solution key");

  // set up the primary variable solution, and its evaluator
  Teuchos::ParameterList pv_sublist = plist_.sublist("primary variable evaluator");
  pv_sublist.set("evaluator name", key_);
  solution_evaluator_ = Teuchos::rcp(new PrimaryVariableFieldEvaluator(pv_sublist));
  S->SetFieldEvaluator(key_, solution_evaluator_);
};


// -----------------------------------------------------------------------------
// Transfer operators -- ONLY COPIES POINTERS
// -----------------------------------------------------------------------------
void PKPhysicalBase::state_to_solution(const Teuchos::RCP<State>& S,
        const Teuchos::RCP<TreeVector>& solution) {
  solution->set_data(S->GetFieldData(key_, name_));
};


// -----------------------------------------------------------------------------
// Transfer operators -- ONLY COPIES POINTERS
// -----------------------------------------------------------------------------
void PKPhysicalBase::solution_to_state(const Teuchos::RCP<TreeVector>& solution,
        const Teuchos::RCP<State>& S) {
  S->SetData(key_, name_, solution->data());
  solution_evaluator_->SetFieldAsChanged();
};


// -----------------------------------------------------------------------------
// Initialization of the PK data.
// -----------------------------------------------------------------------------
void PKPhysicalBase::initialize(const Teuchos::Ptr<State>& S) {
  // initial conditions
  // -- Get the IC function plist.
  if (!plist_.isSublist("initial condition")) {
    std::stringstream messagestream;
    messagestream << name_ << " has no initial condition parameter list.";
    Errors::Message message(messagestream.str());
    Exceptions::amanzi_throw(message);
  }

  // -- Calculate the IC.
  Teuchos::ParameterList ic_plist = plist_.sublist("initial condition");
  Teuchos::RCP<Field> field = S->GetField(key_, name_);
  field->Initialize(ic_plist);

  // -- Push the data into the solution.
  solution_->set_data(field->GetFieldData());
};


} // namespace
