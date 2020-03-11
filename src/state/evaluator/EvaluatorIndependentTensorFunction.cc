/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:

*/

//!

/*
  State

  Author: Ethan Coon

  A field evaluator with no dependencies specified by a function.
*/

#include "EvaluatorIndependentTensorFunction.hh"
#include "UniqueHelpers.hh"

namespace Amanzi {

// ---------------------------------------------------------------------------
// Constructor
// ---------------------------------------------------------------------------
EvaluatorIndependentTensorFunction::EvaluatorIndependentTensorFunction(
  Teuchos::ParameterList& plist)
  : EvaluatorIndependent<TensorVector, TensorVector_Factory>(plist),
    dimension_(-1),
    rank_(-1),
    num_funcs_(-1)
{}

// ---------------------------------------------------------------------------
// Virtual Copy constructor
// ---------------------------------------------------------------------------
Teuchos::RCP<Evaluator>
EvaluatorIndependentTensorFunction::Clone() const
{
  return Teuchos::rcp(new EvaluatorIndependentTensorFunction(*this));
}

// // ---------------------------------------------------------------------------
// // Operator=
// // ---------------------------------------------------------------------------
// Evaluator&
// EvaluatorIndependentTensorFunction::operator=(const Evaluator& other)
// {
//   if (this != &other) {
//     const EvaluatorIndependentTensorFunction* other_p =
//       dynamic_cast<const EvaluatorIndependentTensorFunction*>(&other);
//     AMANZI_ASSERT(other_p != NULL);
//     *this = *other_p;
//   }
//   return *this;
// }

// EvaluatorIndependentTensorFunction&
// EvaluatorIndependentTensorFunction::
// operator=(const EvaluatorIndependentTensorFunction& other)
// {
//   if (this != &other) {
//     AMANZI_ASSERT(my_key_ == other.my_key_);
//     requests_ = other.requests_;
//   }
//   return *this;
// }

void
EvaluatorIndependentTensorFunction::EnsureCompatibility(State& S)
{
  // need only do this once, but AFTER we have a mesh
  auto& f =
    S.Require<TensorVector, TensorVector_Factory>(my_key_, my_tag_, my_key_);
  if (rank_ == -1 && f.map().Mesh().get()) {
    dimension_ = f.dimension();
    AMANZI_ASSERT(dimension_ > 0);
    rank_ = plist_.get<int>("tensor rank");
    f.set_rank(rank_);

    // the map needs to be updated with the correct number of values
    CompositeVectorSpace map_new;
    auto& map_old = f.map();
    map_new.SetMesh(map_old.Mesh());
    num_funcs_ = WhetStone::WHETSTONE_TENSOR_SIZE[dimension_ - 1][rank_ - 1];
    num_funcs_ *= num_funcs_;
    for (auto& name : map_old) {
      map_new.AddComponent(name, map_old.Location(name), num_funcs_);
    }
    f.set_map(map_new);
  }
}

// ---------------------------------------------------------------------------
// Update the value in the state.
// ---------------------------------------------------------------------------
void
EvaluatorIndependentTensorFunction::Update_(State& S)
{
  if (!computed_once_) {
    // Create the function.
    auto& tv = S.Get<TensorVector>(my_key_, my_tag_);
    AMANZI_ASSERT(plist_.isSublist("function"));

    func_ = Functions::createCompositeVectorFunction(plist_.sublist("function"),
            tv.map.Mesh());
  }

  auto& tv = S.GetW<TensorVector>(my_key_, my_tag_, my_key_);
  CompositeVector cv(tv.map.CreateSpace());

  time_ = S.time(my_tag_);
  func_->Compute(time_, cv);
  if (tv.ghosted) cv.ScatterMasterToGhosted();

  // move data into tensor vector
  int j = 0;
  for (auto name : tv.map) {
    auto vec = cv.ViewComponent(name, tv.ghosted);
    Kokkos::parallel_for(
        "TensorFunctionUpdate",
        vec.extent(0),
        KOKKOS_LAMBDA(const int& i) {
          auto Ti = tv.at(j+i);
          Ti(0,0) = vec(i,0);
        });
    j += vec.extent(0);
  }
}

} // namespace Amanzi
