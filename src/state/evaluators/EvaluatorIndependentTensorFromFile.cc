/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon, Bo Gao
*/

/*
  State

*/

#include "EvaluatorIndependentTensorFromFile.hh"
#include "EvaluatorFromFile_Helpers.hh"
#include "Function.hh"
#include "FunctionFactory.hh"

namespace Amanzi {

// ---------------------------------------------------------------------------
// Constructor
// ---------------------------------------------------------------------------
EvaluatorIndependentTensorFromFile::EvaluatorIndependentTensorFromFile(
  Teuchos::ParameterList& plist)
  : EvaluatorIndependent<TensorVector, TensorVector_Factory>(plist),
    filename_(plist.get<std::string>("filename")),
    meshname_(plist.get<std::string>("domain name", "domain")),
    compname_(plist.get<std::string>("component name", "cell")),
    varname_(plist.get<std::string>("variable name")),
    checkpoint_file_(plist.get<bool>("checkpoint file", false)),
    dimension_(-1),
    rank_(-1),
    num_funcs_(-1),
    rescaling_(plist.get<double>("rescaling factor", 1.0))
{
  temporally_variable_ = !plist.get<bool>("constant in time", true);
  if (checkpoint_file_) temporally_variable_ = false;

  if (plist.isParameter("mesh entity")) {
    locname_ = AmanziMesh::createEntityKind(plist.get<std::string>("mesh entity"));
  } else {
    locname_ = AmanziMesh::createEntityKind(compname_);
  }

  if (temporally_variable_ && plist.isSublist("time function")) {
    FunctionFactory fac;
    time_func_ = Teuchos::rcp(fac.Create(plist.sublist("time function")));
  }
}


// ---------------------------------------------------------------------------
// Virtual Copy constructor
// ---------------------------------------------------------------------------
Teuchos::RCP<Evaluator>
EvaluatorIndependentTensorFromFile::Clone() const
{
  return Teuchos::rcp(new EvaluatorIndependentTensorFromFile(*this));
}


// ---------------------------------------------------------------------------
// Operator=
// ---------------------------------------------------------------------------
Evaluator&
EvaluatorIndependentTensorFromFile::operator=(const Evaluator& other)
{
  if (this != &other) {
    const EvaluatorIndependentTensorFromFile* other_p =
      dynamic_cast<const EvaluatorIndependentTensorFromFile*>(&other);
    AMANZI_ASSERT(other_p != NULL);
    *this = *other_p;
  }
  return *this;
}


EvaluatorIndependentTensorFromFile&
EvaluatorIndependentTensorFromFile::operator=(const EvaluatorIndependentTensorFromFile& other)
{
  if (this != &other) {
    AMANZI_ASSERT(my_key_ == other.my_key_);
    requests_ = other.requests_;
  }
  return *this;
}


// ---------------------------------------------------------------------------
// Ensures that the function can provide for the vector's requirements
// ---------------------------------------------------------------------------
void
EvaluatorIndependentTensorFromFile::EnsureCompatibility(State& S)
{
  auto& f = S.Require<TensorVector, TensorVector_Factory>(my_key_, my_tag_, my_key_);

  if (rank_ == -1 && f.map().Mesh().get()) {
    dimension_ = f.dimension();
    AMANZI_ASSERT(dimension_ > 0);

    tensor_type_ = plist_.get<std::string>("tensor type");
    if (tensor_type_ == "scalar") {
      rank_ = 1;
      num_funcs_ = 1;
      ndofs_ = num_funcs_;
    } else if (tensor_type_ == "horizontal and vertical") {
      rank_ = 2;
      num_funcs_ = 2;
      ndofs_ = num_funcs_;
    } else if (tensor_type_ == "diagonal") {
      rank_ = 2;
      num_funcs_ = dimension_;
      ndofs_ = num_funcs_;
    } else if (tensor_type_ == "full symmetric") {
      rank_ = 2;
      num_funcs_ = dimension_ == 2 ? 3 : 6;
      ndofs_ = num_funcs_;
    } else if (tensor_type_ == "full") {
      rank_ = 2;
      num_funcs_ = dimension_ * dimension_;
      ndofs_ = num_funcs_;
    } else {
      Errors::Message msg;
      msg << "EvaluatorIndependentTensorFromFile: invalid parameter, \"tensor type\" = \""
          << tensor_type_
          << "\", must be one of: \"scalar\", \"horizontal and vertical\", \"diagonal\", "
             "\"full symmetric\", or \"full\"";
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

    interpolator_ = Teuchos::rcp(new EvaluatorFromFile_Helpers::FileTimeInterpolator(
      filename_, varname_, compname_, ndofs_, checkpoint_file_));
    interpolator_->setup(map_new, temporally_variable_);
  }
}


// ---------------------------------------------------------------------------
// Update the value in the state
// ---------------------------------------------------------------------------
void
EvaluatorIndependentTensorFromFile::Update_(State& S)
{
  const auto& fac = S.GetRecordSetW(my_key_).GetFactory<TensorVector, TensorVector_Factory>();
  auto& tv = S.GetW<TensorVector>(my_key_, my_tag_, my_key_);

  double t = S.get_time(my_tag_);
  if (time_func_ != Teuchos::null) {
    std::vector<double> point(1, t);
    t = (*time_func_)(point);
  }

  CompositeVector cv_interp(interpolator_->interpolate(t));

  if (rescaling_ != 1.0) cv_interp.Scale(rescaling_);
  if (tv.ghosted) cv_interp.ScatterMasterToGhosted();

  // move data into tensor vector
  int j = 0;
  for (auto name : fac.map()) {
    Epetra_MultiVector& vec = *cv_interp.ViewComponent(name, tv.ghosted);
    EvaluatorFromFile_Helpers::CopyVectorToTensorVector(vec, j, tv);
    j += vec.MyLength();
  }
}


} // namespace Amanzi
