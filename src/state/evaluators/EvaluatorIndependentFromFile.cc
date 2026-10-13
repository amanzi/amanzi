/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

//! An evaluator with no dependencies specified by discrete data in an HDF5
//! file.
/*!

.. todo:
    This needs a test and documentation! --etc

*/

#include "Function.hh"
#include "FunctionFactory.hh"
#include "DataStructuresHelpers.hh"
#include "EvaluatorFromFile_Helpers.hh"
#include "IOFileFormat.hh"

#include "EvaluatorIndependentFromFile.hh"

namespace Amanzi {

const std::string EvaluatorIndependentFromFile::eval_type = "independent variable from file";

// ---------------------------------------------------------------------------
// Constructor
// ---------------------------------------------------------------------------
EvaluatorIndependentFromFile::EvaluatorIndependentFromFile(
  const Teuchos::RCP<Teuchos::ParameterList>& plist)
  : EvaluatorIndependent<CompositeVector, CompositeVectorSpace>(plist),
    filename_(plist->get<std::string>("filename")),
    meshname_(plist->get<std::string>("domain name", "domain")),
    compname_(plist->get<std::string>("component name", "cell")),
    varname_(plist->get<std::string>("variable name")),
    locname_(plist->get<std::string>("mesh entity", "cell")),
    ndofs_(plist->get<int>("number of DoFs", 1)),
    page_size_(plist->get<int>("page size", 2)),
    file_format_(Utils::createIOFileFormat(*plist, filename_)),
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
EvaluatorIndependentFromFile::Clone() const
{
  return Teuchos::rcp(new EvaluatorIndependentFromFile(*this));
}

// ---------------------------------------------------------------------------
// Operator=
// ---------------------------------------------------------------------------
Evaluator&
EvaluatorIndependentFromFile::operator=(const Evaluator& other)
{
  if (this != &other) {
    const EvaluatorIndependentFromFile* other_p =
      dynamic_cast<const EvaluatorIndependentFromFile*>(&other);
    AMANZI_ASSERT(other_p != NULL);
    *this = *other_p;
  }
  return *this;
}

EvaluatorIndependentFromFile&
EvaluatorIndependentFromFile::operator=(const EvaluatorIndependentFromFile& other)
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
EvaluatorIndependentFromFile::EnsureCompatibility(State& S)
{
  EvaluatorIndependent::EnsureCompatibility(S);

  if (interpolator_ == Teuchos::null) {
    auto mesh = S.GetMesh(meshname_);

    // requirements on vector data
    AmanziMesh::Entity_kind location = AmanziMesh::createEntityKind(locname_);

    auto& space = S.Require<CompositeVector, CompositeVectorSpace>(my_key_, my_tag_, my_key_);
    space.SetMesh(mesh)->AddComponent(compname_, location, ndofs_);

    interpolator_ =
      Teuchos::rcp(new EvaluatorFromFile_Helpers::TimeSeriesInterpolatorCV(filename_,
                                                                           varname_,
                                                                           space.CreateSpace(),
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
EvaluatorIndependentFromFile::Update_(State& S)
{
  CompositeVector& cv = S.GetW<CompositeVector>(my_key_, my_tag_, my_key_);

  double t = S.get_time(my_tag_);
  if (time_func_ != Teuchos::null) {
    Kokkos::View<double**, Kokkos::HostSpace> point("time", 1, 1);
    point(0, 0) = t;
    t = (*time_func_)(point);
  }

  interpolator_->interpolate(t, cv);

  if (locname_ == "cell" && (cv.hasComponent("boundary_face") || cv.hasComponent("face")))
    DeriveFaceValuesFromCellValues(cv);
}

} // namespace Amanzi
