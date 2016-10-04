/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------
   ATS

   License: see $ATS_DIR/COPYRIGHT
   Author: Ethan Coon

   Default base with default implementations of methods for a physical PK.
   ------------------------------------------------------------------------- */
#include "pk_physical_base.hh"

namespace Amanzi {


<<<<<<< HEAD
  PKPhysicalBase::PKPhysicalBase(Teuchos::Ptr<State> S, const Teuchos::RCP<Teuchos::ParameterList>& plist,
        Teuchos::ParameterList& FElist,
        const Teuchos::RCP<TreeVector>& solution) :
    PKDefaultBase(S, plist,FElist,solution) {
=======
PKPhysicalBase::PKPhysicalBase(const Teuchos::RCP<Teuchos::ParameterList>& plist,
        Teuchos::ParameterList& FElist,
        const Teuchos::RCP<TreeVector>& solution) :
    PKDefaultBase(plist,FElist,solution) {
>>>>>>> 3712d1ddeb1cfe9f074d84ba39b930e7f970357e

  // domain -- default is the entire mesh, no prefix
  if (domain_.empty()) {
   domain_ = plist_->get<std::string>("domain name", std::string("domain"));
  }
  
  if (key_.empty()) {
    key_ = plist_->get<std::string>("primary variable");
  }

  // set up the primary variable solution, and its evaluator
  Teuchos::ParameterList& pv_sublist = FElist.sublist(key_);
  pv_sublist.set("evaluator name", key_);
  pv_sublist.set("field evaluator type", "primary variable");

  // primary variable max change
  max_valid_change_ = plist_->get<double>("max valid change", -1.0);

}

// -----------------------------------------------------------------------------
// Construction of data.
// -----------------------------------------------------------------------------
void PKPhysicalBase::setup(const Teuchos::Ptr<State>& S) {
  PKDefaultBase::setup(S);

  // get the mesh
  mesh_ = S->GetMesh(domain_);

<<<<<<< HEAD

=======
>>>>>>> 3712d1ddeb1cfe9f074d84ba39b930e7f970357e
  // set up the debugger
  db_ = Teuchos::rcp(new Debugger(mesh_, name_, *plist_));

  // require primary variable evaluator
  S->RequireFieldEvaluator(key_);
  Teuchos::RCP<FieldEvaluator> fm = S->GetFieldEvaluator(key_);
  solution_evaluator_ = Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>(fm);
  ASSERT(solution_evaluator_ != Teuchos::null);
};


// -----------------------------------------------------------------------------
// Transfer operators -- ONLY COPIES POINTERS
// -----------------------------------------------------------------------------
void PKPhysicalBase::state_to_solution(const Teuchos::RCP<State>& S,
        TreeVector& solution) {
  solution.SetData(S->GetFieldData(key_, name_));
};


// -----------------------------------------------------------------------------
// Transfer operators -- ONLY COPIES POINTERS
// -----------------------------------------------------------------------------
void PKPhysicalBase::solution_to_state(TreeVector& solution,
        const Teuchos::RCP<State>& S) {
  ASSERT(solution.Data() == S->GetFieldData(key_));
  //  S->SetData(key_, name_, solution->Data());
  //  solution_evaluator_->SetFieldAsChanged();
};


// -----------------------------------------------------------------------------
// Experimental approach -- we must pull out S_next_'s solution_evaluator_ to
// stay current for ChangedSolution()
// -----------------------------------------------------------------------------
void PKPhysicalBase::set_states(const Teuchos::RCP<const State>& S,
        const Teuchos::RCP<State>& S_inter,
        const Teuchos::RCP<State>& S_next) {
  PKDefaultBase::set_states(S, S_inter, S_next);

  // Get the FE and mark it as changed.
  // Note that this is necessary because we need this to point at the
  // FE in S_next_, not the one which we created in S_.
  Teuchos::RCP<FieldEvaluator> fm = S_next->GetFieldEvaluator(key_);
#if ENABLE_DBC
  solution_evaluator_ = Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>(fm);
  ASSERT(solution_evaluator_ != Teuchos::null);
#else
  solution_evaluator_ = Teuchos::rcp_static_cast<PrimaryVariableFieldEvaluator>(fm);
#endif

  solution_evaluator_->SetFieldAsChanged(S_next_.ptr());
};


// -----------------------------------------------------------------------------
// Ensures the step size is smaller than max_valid_change
// -----------------------------------------------------------------------------
bool PKPhysicalBase::valid_step() {
  Teuchos::OSTab tab = vo_->getOSTab();
  if (vo_->os_OK(Teuchos::VERB_EXTREME))
    *vo_->os() << "Validating time step." << std::endl;

  if (max_valid_change_ > 0.0) {
    const CompositeVector& var_new = *S_next_->GetFieldData(key_);
    const CompositeVector& var_old = *S_->GetFieldData(key_);
    CompositeVector dvar(var_new);
    dvar.Update(-1., var_old, 1.);
    double change = 0.;
    dvar.NormInf(&change);
    if (change > max_valid_change_) {
      if (vo_->os_OK(Teuchos::VERB_MEDIUM))
        *vo_->os() << "Invalid time step, max primary variable change="
                   << change << " > limit=" << max_valid_change_ << std::endl;
      return false;
    }
  }
  return true;
}


// -----------------------------------------------------------------------------
// Initialization of the PK data.
// -----------------------------------------------------------------------------
void PKPhysicalBase::initialize(const Teuchos::Ptr<State>& S) {
  Teuchos::RCP<Field> field = S->GetField(key_, name_);

  if (!field->initialized()) {
    // initial conditions
    // -- Get the IC function plist.
    if (!plist_->isSublist("initial condition")) {
      std::stringstream messagestream;
      messagestream << name_ << " has no initial condition parameter list.";
      Errors::Message message(messagestream.str());
      Exceptions::amanzi_throw(message);
    }

    // -- Calculate the IC.
    Teuchos::ParameterList ic_plist = plist_->sublist("initial condition");
    field->Initialize(ic_plist);
<<<<<<< HEAD
=======

>>>>>>> 3712d1ddeb1cfe9f074d84ba39b930e7f970357e
    // -- Update faces from cells if needed.
    // if (ic_plist.get<bool>("initialize faces from cells", false)) {
    //   DeriveFaceValuesFromCellValues_(field->GetFieldData().ptr());
    // }

    // communicate just to make sure values are initialized for valgrind's sake
    field->GetFieldData()->ScatterMasterToGhosted();
    solution_evaluator_->SetFieldAsChanged(S);
  }

  // -- Push the data into the solution.
  solution_->SetData(field->GetFieldData());
};


// -----------------------------------------------------------------------------
// Interpolate pressure ICs on cells to ICs for lambda (faces).
// -----------------------------------------------------------------------------
void PKPhysicalBase::DeriveFaceValuesFromCellValues_(const Teuchos::Ptr<CompositeVector>& cv) {
  cv->ScatterMasterToGhosted("cell");
  Teuchos::Ptr<const CompositeVector> cv_const(cv);
  const Epetra_MultiVector& cv_c = *cv_const->ViewComponent("cell",true);
  Epetra_MultiVector& cv_f = *cv->ViewComponent("face",false);

  int f_owned = cv_f.MyLength();
  for (int f=0; f!=f_owned; ++f) {
    AmanziMesh::Entity_ID_List cells;
    cv->Mesh()->face_get_cells(f, AmanziMesh::USED, &cells);
    int ncells = cells.size();

    double face_value = 0.0;
    for (int n=0; n!=ncells; ++n) {
      face_value += cv_c[0][cells[n]];
    }
    cv_f[0][f] = face_value / ncells;
  }
};


} // namespace
