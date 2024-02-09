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

#include "CompositeVectorFunction.hh"
#include "EvaluatorIndependentTensorFunction.hh"

namespace Amanzi {

// ---------------------------------------------------------------------------
// Constructor
// ---------------------------------------------------------------------------
EvaluatorIndependentTensorFunction::EvaluatorIndependentTensorFunction(
  const Teuchos::RCP<Teuchos::ParameterList>& plist)
  : EvaluatorIndependent<TensorVector, TensorVector_Factory>(plist),
    dimension_(-1),
    rank_(-1),
    num_funcs_(-1),
    rescaling_(plist->get<double>("rescaling factor", 1.0))
{}

// ---------------------------------------------------------------------------
// Virtual Copy constructor
// ---------------------------------------------------------------------------
Teuchos::RCP<Evaluator>
EvaluatorIndependentTensorFunction::Clone() const
{
  return Teuchos::rcp(new EvaluatorIndependentTensorFunction(*this));
}


void
EvaluatorIndependentTensorFunction::EnsureCompatibility(State& S)
{
  // need only do this once, but AFTER we have a mesh
  auto& f =
    S.Require<TensorVector, TensorVector_Factory>(my_key_, my_tag_, my_key_);
  if (rank_ == -1 && f.getMap().getMesh().get()) {
    dimension_ = f.getDimension();
    AMANZI_ASSERT(dimension_ > 0);
    rank_ = plist_->get<int>("tensor rank");
    if (rank_ < 1 || rank_ > 4 || rank_ == 3) {
      Errors::Message msg;
      msg << "EvaluatorIndependentTensorFunction: Invalid \"tensor rank\" of " << rank_ << " requested, valid are 1 (scalar), 2 (square), and 4.";
      Exceptions::amanzi_throw(msg);
    }
    f.setRank(rank_);

    // the map needs to be updated with the correct number of values
    CompositeVectorSpace map_new;
    auto& map_old = f.getMap();
    map_new.SetMesh(map_old.getMesh());

    // this is not correct -- should be this for a generic, nonsymmetric
    // tensor.  In practice we will want to supply options to have, e.g. rank
    // but also diagonal (but not isotropic), symmetric, horizontal + vertical,
    // etc
    int tensor_size = WhetStone::WHETSTONE_TENSOR_SIZE[dimension_ - 1][rank_ - 1];
    num_funcs_ = tensor_size * tensor_size;
    for (auto& name : map_old) {
      map_new.AddComponent(name, map_old.getLocation(name), num_funcs_);
    }
    f.setMap(map_new);
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
    AMANZI_ASSERT(plist_->isSublist("function"));

    func_ = Teuchos::rcp(new Functions::CompositeVectorFunction(plist_->sublist("function"), tv.map.getMesh()));
  }

  auto& tv = S.GetW<TensorVector>(my_key_, my_tag_, my_key_);
  CompositeVector cv(tv.map.CreateSpace());

  time_ = S.get_time(my_tag_);
  func_->Compute(time_, cv);
  if (rescaling_ != 1.0) cv.scale(rescaling_);
  if (tv.ghosted) cv.scatterMasterToGhosted();

  // move data into tensor vector
  int j = 0;
  for (auto name : tv.map) {
    auto vec = cv.viewComponent(name, tv.ghosted);
    Impl::assignViewToTensorVectorDiag(vec, j, tv);
    j += vec.extent(0);
  }
}

} // namespace Amanzi
