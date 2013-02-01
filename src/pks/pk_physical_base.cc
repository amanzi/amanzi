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

  // If they have not been set, pull the domain name and primary variable key
  // from the parameter list.

  // domain
  if (domain_ == std::string("")) {
    domain_ = plist_.get<std::string>("domain name", std::string("domain"));
  }
  if (key_ == std::string("")) {
    key_ = plist_.get<std::string>("primary variable key");
  }

  // derive the prefix
  if (domain_ == std::string("domain")) {
    domain_prefix_ = std::string("");
  } else {
    domain_prefix_ = domain_ + std::string("_");
  }

  // get the mesh
  mesh_ = S->GetMesh(domain_);

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
  //  solution_evaluator_->SetFieldAsChanged();
};


// -----------------------------------------------------------------------------
// Initialization of the PK data.
// -----------------------------------------------------------------------------
void PKPhysicalBase::initialize(const Teuchos::Ptr<State>& S) {
  Teuchos::RCP<Field> field = S->GetField(key_, name_);

  if (!field->initialized()) {
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
    field->Initialize(ic_plist);

    // -- Update faces from cells if needed.
    if (ic_plist.get<bool>("initialize faces from cells", false)) {
      DeriveFaceValuesFromCellValues_(field->GetFieldData().ptr());
    }
  }

  // -- Push the data into the solution.
  solution_->set_data(field->GetFieldData());
};


// -----------------------------------------------------------------------------
// Interpolate pressure ICs on cells to ICs for lambda (faces).
// -----------------------------------------------------------------------------
void PKPhysicalBase::DeriveFaceValuesFromCellValues_(const Teuchos::Ptr<CompositeVector>& cv) {
  AmanziMesh::Entity_ID_List cells;
  cv->ScatterMasterToGhosted("cell");

  int f_owned = cv->size("face");
  for (int f=0; f!=f_owned; ++f) {
    cells.clear();
    cv->mesh()->face_get_cells(f, AmanziMesh::USED, &cells);
    int ncells = cells.size();

    double face_value = 0.0;
    for (int n=0; n!=ncells; ++n) {
      face_value += (*cv)("cell",cells[n]);
    }
    (*cv)("face",f) = face_value / ncells;
  }
};


} // namespace
