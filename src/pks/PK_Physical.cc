/* -------------------------------------------------------------------------
  Process Kernels

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Daniil Svyatsky
*/

#include "PK_Physical.hh"

#include "EvaluatorIndependentFunction.hh"
#include "State.hh"
#include "TreeVector.hh"

namespace Amanzi {

// -----------------------------------------------------------------------------
// Transfer operators -- copies ONLY pointers
// -----------------------------------------------------------------------------
void
PK_Physical::State_to_Solution(const Tag& tag, TreeVector& solution)
{
  solution.SetData(S_->GetPtrW<CompositeVector>(key_, tag, name()));
}


void
PK_Physical::Solution_to_State(const TreeVector& solution, const Tag& tag)
{
  AMANZI_ASSERT(solution.Data() == S_->GetPtr<CompositeVector>(key_, tag));
}


// -----------------------------------------------------------------------------
// Helper method to add a primary variable evaluator
// -----------------------------------------------------------------------------
void
PK_Physical::AddDefaultPrimaryEvaluator_(const Key& key, const Tag& tag)
{
  AMANZI_ASSERT(S_ != Teuchos::null);
  Teuchos::ParameterList elist = S_->GetEvaluatorList(key);
  elist.set<std::string>("tag", tag.get());
  elist.setName(key);
  auto eval = Teuchos::rcp(new EvaluatorPrimary<CompositeVector, CompositeVectorSpace>(elist));
  S_->SetEvaluator(key, tag, eval);
}


// -----------------------------------------------------------------------------
// Helper method to add an independent variable evaluator
// -----------------------------------------------------------------------------
void
PK_Physical::AddDefaultIndependentEvaluator_(const Key& key, const Tag& tag, double val)
{
  Teuchos::ParameterList elist(key);
  elist.set<std::string>("evaluator type", "independent variable")
    .sublist("function")
    .sublist("ALL")
    .set<std::string>("region", "All")
    .set<std::string>("component", "*")
    .sublist("function")
    .sublist("function-constant")
    .set<double>("value", val);

  auto eval = Teuchos::rcp(new EvaluatorIndependentFunction(elist));
  S_->SetEvaluator(key, tag, eval);
}


// -----------------------------------------------------------------------------
// Helper method to initialize a CV field
// -----------------------------------------------------------------------------
void
InitializeCVField(const Teuchos::RCP<State>& S,
                  const VerboseObject& vo,
                  const Key& key,
                  const Tag& tag,
                  const Key& passwd,
                  double default_val)
{
  Teuchos::OSTab tab = vo.getOSTab();

  if (S->HasRecord(key, tag)) {
    if (S->GetRecord(key, tag).owner() == passwd) {
      if (!S->GetRecord(key, tag).initialized()) {
        S->GetW<CompositeVector>(key, tag, passwd).PutScalar(default_val);
        S->GetRecordW(key, tag, passwd).set_initialized();

        if (vo.os_OK(Teuchos::VERB_MEDIUM))
          *vo.os() << "initialized \"" << key << "\" to value " << default_val << std::endl;
      }
    }
  }
}

} // namespace Amanzi
