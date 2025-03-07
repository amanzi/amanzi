/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

#include "EvaluatorTemporalInterpolation.hh"

namespace Amanzi {

// ---------------------------------------------------------------------------
// Constructor
// ---------------------------------------------------------------------------
EvaluatorTemporalInterpolation::EvaluatorTemporalInterpolation(Teuchos::ParameterList& plist)
  : EvaluatorSecondary(plist)
{
  Key my_key = my_keys_.front().first;

  Tag current_tag(plist.get<std::string>("current tag"));
  current_ = KeyTag{ my_key, current_tag };
  dependencies_.insert(current_);

  Tag next_tag(plist.get<std::string>("next tag"));
  next_ = KeyTag{ my_key, next_tag };
  dependencies_.insert(next_);

  // also depend upon time at all three tags
  dependencies_.insert(KeyTag{"time", my_keys_.front().second});
  dependencies_.insert(KeyTag{"time", current_tag});
  dependencies_.insert(KeyTag{"time", next_tag});
}

// ---------------------------------------------------------------------------
// Copy constructor
// ---------------------------------------------------------------------------
Teuchos::RCP<Evaluator>
EvaluatorTemporalInterpolation::Clone() const //override
{
  return Teuchos::rcp(new EvaluatorTemporalInterpolation(*this));
}


void
EvaluatorTemporalInterpolation::EnsureCompatibility(State& S)
{
  // claim ownership
  auto& fac = S.Require<CompositeVector, CompositeVectorSpace>(my_keys_.front().first,
          my_keys_.front().second, my_keys_.front().first);

  EnsureCompatibility_Flags_(S);

  // structure the same across all vectors
  for (const auto& dep : dependencies_) {
    if (dep.first == "time") {
      S.Require<double>(dep.first, dep.second);
    } else {
      S.Require<CompositeVector, CompositeVectorSpace>(dep.first, dep.second).Update(fac);
    }
  }
}


void
EvaluatorTemporalInterpolation::Update_(State& S)
{
  double t_current = S.get_time(current_.second);
  double t_next = S.get_time(next_.second);
  double t_this = S.get_time(my_keys_.front().second);

  CompositeVector& result = S.GetW<CompositeVector>(my_keys_.front().first, my_keys_.front().second, my_keys_.front().first);
  const CompositeVector& var_current = S.Get<CompositeVector>(current_.first, current_.second);
  const CompositeVector& var_next = S.Get<CompositeVector>(next_.first, next_.second);
  AMANZI_ASSERT(t_this >= t_current);
  AMANZI_ASSERT(t_this <= t_next);

  double a(1); // if dt is 0, default to next value
  if (t_next - t_current > 0.) {
    a = (t_this - t_current) / (t_next - t_current);
  }
  result.Update(1-a, var_current, a, var_next, 0);
}

// void
// EvaluatorTemporalInterpolation::EvaluatePartialDerivative_(
//         const State& S,
//         const Key& wrt_key,
//         const Tag& wrt_tag,
//         const std::vector<CompositeVector*>& result)
// {
//   AMANZI_ASSERT(dependencies_.size() == 2);
//   double t_current = S.get_time(current_.second);
//   double t_next = S.get_time(next_.second);
//   double t_this = S.get_time(my_keys_.front().second);
//   AMANZI_ASSERT(t_this >= t_current);
//   AMANZI_ASSERT(t_this <= t_next);

//   double a(1.);
//   if (t_next - t_current > 0.) {
//     a = (t_this - t_current) / (t_next - t_current);
//   }
//   if (wrt_tag == current_.second) {
//     result[0]->PutScalar(1-a);
//   } else if (wrt_tag == next_.second) {
//     result[0]->PutScalar(a);
//   } else {
//     AMANZI_ASSERT(false);
//   }
// }

} // namespace Amanzi
