/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/*
  State

  A field evaluator with no dependencies specified by a function.
*/

#include <memory>

#include "EvaluatorIndependentFunction.hh"
#include "CompositeVectorFunctionFactory.hh"

namespace Amanzi {

// ---------------------------------------------------------------------------
// Constructor
// ---------------------------------------------------------------------------
EvaluatorIndependentFunction::EvaluatorIndependentFunction(Teuchos::ParameterList& plist)
  : EvaluatorIndependent<CompositeVector, CompositeVectorSpace>(plist),
    dot_with_normal_(plist.get<bool>("dot with normal", false)),
    spatial_dist_method_(plist.get<std::string>("spatial distribution method", "none")) {};


// ---------------------------------------------------------------------------
// Virtual Copy constructor
// ---------------------------------------------------------------------------
Teuchos::RCP<Evaluator>
EvaluatorIndependentFunction::Clone() const
{
  return Teuchos::rcp(new EvaluatorIndependentFunction(*this));
}


// ---------------------------------------------------------------------------
// Operator=
// ---------------------------------------------------------------------------
Evaluator&
EvaluatorIndependentFunction::operator=(const Evaluator& other)
{
  if (this != &other) {
    const EvaluatorIndependentFunction* other_p =
      dynamic_cast<const EvaluatorIndependentFunction*>(&other);
    AMANZI_ASSERT(other_p != NULL);
    *this = *other_p;
  }
  return *this;
}


EvaluatorIndependentFunction&
EvaluatorIndependentFunction::operator=(const EvaluatorIndependentFunction& other)
{
  if (this != &other) {
    AMANZI_ASSERT(my_key_ == other.my_key_);
    requests_ = other.requests_;
  }
  return *this;
}


void
EvaluatorIndependentFunction::EnsureCompatibility(State& S)
{
  EvaluatorIndependent<CompositeVector, CompositeVectorSpace>::EnsureCompatibility(S);

  auto& fac = S.Require<CompositeVector, CompositeVectorSpace>(my_key_, my_tag_);
  if (fac.Mesh() != Teuchos::null) {
    std::vector<std::string> complist;
    func_ = Functions::CreateCompositeVectorFunction(
      plist_.sublist("function"), fac, complist, dot_with_normal_, spatial_dist_method_);
  }
}


// ---------------------------------------------------------------------------
// Update the value in the state.
// ---------------------------------------------------------------------------
void
EvaluatorIndependentFunction::Update_(State& S)
{
  // NOTE: EvaluatorIndependentFunctions own their own data.
  auto cv = S.GetPtrW<CompositeVector>(my_key_, my_tag_, my_key_);
  time_ = S.get_time(my_tag_);
  func_->Compute(time_, cv.ptr());
}

} // namespace Amanzi
