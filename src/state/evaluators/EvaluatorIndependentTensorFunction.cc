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

#include "CompositeVectorFunctionFactory.hh"
#include "EvaluatorIndependentTensorFunction.hh"
#include "EvaluatorFromFile_Helpers.hh"

namespace Amanzi {


// ---------------------------------------------------------------------------
// Constructor
// ---------------------------------------------------------------------------
EvaluatorIndependentTensorFunction::EvaluatorIndependentTensorFunction(
  Teuchos::ParameterList& plist)
  : EvaluatorIndependent<TensorVector, TensorVector_Factory>(plist),
    dimension_(-1),
    rank_(-1),
    num_funcs_(-1),
    rescaling_(plist.get<double>("rescaling factor", 1.0))
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
  auto& f = S.Require<TensorVector, TensorVector_Factory>(my_key_, my_tag_, my_key_);
  if (rank_ == -1 && f.map().Mesh().get()) {
    dimension_ = f.dimension();
    AMANZI_ASSERT(dimension_ > 0);

    tensor_type_ = plist_.get<std::string>("tensor type");
    if (tensor_type_ == "scalar") {
      rank_ = 1;
      num_funcs_ = 1;
    } else if (tensor_type_ == "horizontal and vertical") {
      rank_ = 2;
      num_funcs_ = 2;
    } else if (tensor_type_ == "diagonal") {
      rank_ = 2;
      num_funcs_ = dimension_;
    } else if (tensor_type_ == "full symmetric") {
      rank_ = 2;
      num_funcs_ = dimension_ == 2 ? 3 : 6;
    } else if (tensor_type_ == "full") {
      rank_ = 2;
      num_funcs_ = dimension_ * dimension_;
    } else {
      Errors::Message msg;
      msg << "EvaluatorIndependentTensorFunction: invalid paramter, \"tensor type\" = \""
          << tensor_type_
          << "\", must be one of: \"scalar\", \"horizontal and vertical\", \"diagonal\", \"full "
             "symmetric\""
          << ", or \"full\"";
      Exceptions::amanzi_throw(msg);
    }

    f.set_rank(rank_);

    // the map needs to be updated with the correct number of values
    CompositeVectorSpace map_new;
    auto& map_old = f.map();
    map_new.SetMesh(map_old.Mesh());

    for (auto& name : map_old) {
      map_new.AddComponent(name, map_old.Location(name), num_funcs_);
    }
    f.set_map(map_new);

    std::vector<std::string> component_names;
    func_ = Functions::CreateCompositeVectorFunction(
      plist_.sublist("function"), map_new, component_names);
  }
}

// ---------------------------------------------------------------------------
// Update the value in the state.
// ---------------------------------------------------------------------------
void
EvaluatorIndependentTensorFunction::Update_(State& S)
{
  const auto& fac = S.GetRecordSetW(my_key_).GetFactory<TensorVector, TensorVector_Factory>();
  // if (!computed_once_) {
  //   // Create the function.
  //   AMANZI_ASSERT(plist_.isSublist("function"));

  //   std::vector<std::string> component_names;
  //   func_ = Functions::CreateCompositeVectorFunction(plist_.sublist("function"), fac.map(), component_names);
  // }

  auto& tv = S.GetW<TensorVector>(my_key_, my_tag_, my_key_);
  auto cv = Teuchos::rcp(new CompositeVector(fac.map()));

  time_ = S.get_time(my_tag_);
  func_->Compute(time_, cv.ptr());

  if (rescaling_ != 1.0) cv->Scale(rescaling_);
  if (tv.ghosted) cv->ScatterMasterToGhosted();

  // move data into tensor vector
  int j = 0;
  for (auto name : fac.map()) {
    Epetra_MultiVector& vec = *cv->ViewComponent(name, tv.ghosted);
    EvaluatorFromFile_Helpers::CopyVectorToTensorVector(vec, j, tv);
    j += vec.MyLength();
  }
}


} // namespace Amanzi
