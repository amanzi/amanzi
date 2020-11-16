/* -------------------------------------------------------------------------
  Process Kernels

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Daniil Svyatsky
*/

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
}


// -----------------------------------------------------------------------------
// Transfer operators -- ONLY COPIES POINTERS
// -----------------------------------------------------------------------------
void PK_Physical::Solution_to_State(TreeVector& solution,
                                    const Teuchos::RCP<State>& S) {
  AMANZI_ASSERT(solution.Data() == S->GetFieldData(key_));
  //  S->SetData(key_, name_, solution->Data());
  //  solution_evaluator_->SetFieldAsChanged();
}


void PK_Physical::Solution_to_State(const TreeVector& solution,
                                    const Teuchos::RCP<State>& S) {
  AMANZI_ASSERT(solution.Data() == S->GetFieldData(key_));
  //  S->SetData(key_, name_, solution->Data());
  //  solution_evaluator_->SetFieldAsChanged();
}
  

// -----------------------------------------------------------------------------
// Experimental approach -- we must pull out S_next_'s solution_evaluator_ to
// stay current for ChangedSolution()
// -----------------------------------------------------------------------------
void PK_Physical::set_states(const Teuchos::RCP<State>& S,
                             const Teuchos::RCP<State>& S_inter,
                             const Teuchos::RCP<State>& S_next) {
  S_ = S;
  S_inter_ = S_inter;
  S_next_ = S_next;
    
  // Get the FE and mark it as changed.
  // Note that this is necessary because we need this to point at the
  // FE in S_next_, not the one which we created in S_.
  Teuchos::RCP<FieldEvaluator> fm = S_next->GetFieldEvaluator(key_);
#if ENABLE_DBC
  solution_evaluator_ = Teuchos::rcp_dynamic_cast<PrimaryVariableFieldEvaluator>(fm);
  AMANZI_ASSERT(solution_evaluator_ != Teuchos::null);
#else
  solution_evaluator_ = Teuchos::rcp_static_cast<PrimaryVariableFieldEvaluator>(fm);
#endif

  solution_evaluator_->SetFieldAsChanged(S_next_.ptr());
}


// -----------------------------------------------------------------------------
// Populate state field
// -----------------------------------------------------------------------------
void PK_Physical::InitializeField_(const Teuchos::Ptr<State>& S, 
                                  const std::string& passwd,
                                  std::string fieldname, double default_val)
{
  Teuchos::OSTab tab = vo_->getOSTab();

  if (S->HasField(fieldname)) {
    if (S->GetField(fieldname)->owner() == passwd) {
      if (!S->GetField(fieldname)->initialized()) {
        S->GetFieldData(fieldname, passwd)->PutScalar(default_val);
        S->GetField(fieldname, passwd)->set_initialized();

        if (vo_->os_OK(Teuchos::VERB_MEDIUM))
           *vo_->os() << "initialized \"" << fieldname << "\" to value " << default_val << std::endl;
      }
    }
  }
}


void PK_Physical::AddDefaultPrimaryEvaluator_(const Key& key)
{
  Teuchos::ParameterList elist;
  elist.set<std::string>("evaluator name", key);
  auto eval = Teuchos::rcp(new PrimaryVariableFieldEvaluator(elist));
  AMANZI_ASSERT(S_ != Teuchos::null);
  S_->SetFieldEvaluator(key, eval);
}

} // namespace Amanzi

