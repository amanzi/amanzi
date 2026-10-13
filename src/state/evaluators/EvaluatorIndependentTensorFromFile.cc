/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

//! An evaluator for a WhetStone::Tensor field, specified by discrete data in
//! an HDF5 file.
/*!

.. todo:
    This needs a test! --etc

*/

#include "Function.hh"
#include "FunctionFactory.hh"
#include "EvaluatorFromFile_Helpers.hh"
#include "EvaluatorTensorFromFile_Helpers.hh"
#include "IOFileFormat.hh"

#include "EvaluatorIndependentTensorFromFile.hh"

namespace Amanzi {

const std::string EvaluatorIndependentTensorFromFile::eval_type =
  "independent variable tensor from file";

// ---------------------------------------------------------------------------
// Constructor
// ---------------------------------------------------------------------------
EvaluatorIndependentTensorFromFile::EvaluatorIndependentTensorFromFile(
  const Teuchos::RCP<Teuchos::ParameterList>& plist)
  : EvaluatorIndependent<TensorVector, TensorVector_Factory>(plist),
    filename_(plist->get<std::string>("filename")),
    meshname_(plist->get<std::string>("domain name", "domain")),
    compname_(plist->get<std::string>("component name", "cell")),
    varname_(plist->get<std::string>("variable name")),
    locname_(plist->get<std::string>("mesh entity", "cell")),
    tensor_type_(plist->get<std::string>("tensor type")),
    dimension_(-1),
    rank_(-1),
    num_funcs_(-1),
    page_size_(plist->get<int>("page size", 2)),
    file_format_(Utils::createIOFileFormat(*plist, filename_)),
    rescaling_(plist->get<double>("rescaling factor", 1.0)),
    checkpoint_file_(plist->get<bool>("checkpoint file", false))
{
  if (checkpoint_file_) temporally_variable_ = false;

  if (temporally_variable_ && plist->isSublist("time function")) {
    FunctionFactory fac;
    time_func_ = Teuchos::rcp(fac.Create(plist->sublist("time function")));
  }
}

// ---------------------------------------------------------------------------
// Virtual Copy constructor
// ---------------------------------------------------------------------------
Teuchos::RCP<Evaluator>
EvaluatorIndependentTensorFromFile::Clone() const
{ return Teuchos::rcp(new EvaluatorIndependentTensorFromFile(*this)); }

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
// Ensures that the function can provide for the vector's requirements.
// ---------------------------------------------------------------------------
void
EvaluatorIndependentTensorFromFile::EnsureCompatibility(State& S)
{
  // need only do this once, but AFTER we have a mesh
  auto& f = S.Require<TensorVector, TensorVector_Factory>(my_key_, my_tag_, my_key_);
  if (rank_ == -1 && f.getMap().getMesh().get()) {
    dimension_ = f.getDimension();
    AMANZI_ASSERT(dimension_ > 0);

    EvaluatorTensorFromFile_Helpers::parseTensorType(
      tensor_type_, dimension_, "EvaluatorIndependentTensorFromFile", rank_, num_funcs_);
    f.setRank(rank_);

    // the map needs to be updated with the correct number of values (DoFs)
    CompositeVectorSpace map_new;
    auto& map_old = f.getMap();
    map_new.SetMesh(map_old.getMesh());

    AmanziMesh::Entity_kind location = AmanziMesh::createEntityKind(locname_);
    map_new.AddComponent(compname_, location, num_funcs_);
    f.setMap(map_new);

    interpolator_ =
      Teuchos::rcp(new EvaluatorFromFile_Helpers::TimeSeriesInterpolatorCV(filename_,
                                                                           varname_,
                                                                           map_new.CreateSpace(),
                                                                           temporally_variable_,
                                                                           page_size_,
                                                                           file_format_,
                                                                           checkpoint_file_));
  }
}

// ---------------------------------------------------------------------------
// Update the value in the state.
// ---------------------------------------------------------------------------
void
EvaluatorIndependentTensorFromFile::Update_(State& S)
{
  const auto& fac = S.GetFactoryW<TensorVector, TensorVector_Factory>(my_key_);
  auto& tv = S.GetW<TensorVector>(my_key_, my_tag_, my_key_);

  double t = S.get_time(my_tag_);
  if (time_func_ != Teuchos::null) {
    Kokkos::View<double**, Kokkos::HostSpace> point("time", 1, 1);
    point(0, 0) = t;
    t = (*time_func_)(point);
  }

  CompositeVector cv(interpolator_->getSpace());
  interpolator_->interpolate(t, cv);

  if (rescaling_ != 1.0) cv.scale(rescaling_);
  if (tv.ghosted) cv.scatterMasterToGhosted();

  // move data into tensor vector, on the host, then sync to device
  int j = 0;
  for (auto name : fac->getMap()) {
    auto vec = cv.viewComponent<MemSpace_kind::HOST>(name, tv.ghosted);
    EvaluatorTensorFromFile_Helpers::copyVectorToTensorVector(vec, j, dimension_, tv);
    j += vec.extent(0);
  }
  tv.data.update_entries_device();
}

} // namespace Amanzi
