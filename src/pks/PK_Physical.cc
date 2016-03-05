/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/* -------------------------------------------------------------------------


Author: Daniil Svyatsky


------------------------------------------------------------------------- */

#include "PK_Physical.hh"
#include "State.hh"
#include "TreeVector.hh"

namespace Amanzi {

  // -----------------------------------------------------------------------------
  // Transfer operators -- ONLY COPIES POINTERS
  // -----------------------------------------------------------------------------
  void PK_Physical::State_to_Solution(const Teuchos::RCP<State>& S,
                                         TreeVector& solution) {
    solution.SetData(S->GetFieldData(key_, name_));
  };


  // -----------------------------------------------------------------------------
  // Transfer operators -- ONLY COPIES POINTERS
  // -----------------------------------------------------------------------------
  void PK_Physical::Solution_to_State(TreeVector& solution,
                                         const Teuchos::RCP<State>& S) {
    ASSERT(solution.Data() == S->GetFieldData(key_));
    //  S->SetData(key_, name_, solution->Data());
    //  solution_evaluator_->SetFieldAsChanged();
  };


  // -----------------------------------------------------------------------------
  // Experimental approach -- we must pull out S_next_'s solution_evaluator_ to
  // stay current for ChangedSolution()
  // -----------------------------------------------------------------------------
  void PK_Physical::set_states(const Teuchos::RCP<const State>& S,
                                  const Teuchos::RCP<State>& S_inter,
                                  const Teuchos::RCP<State>& S_next) {
    PK_Default::set_states(S, S_inter, S_next);

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


}
