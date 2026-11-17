/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/*
  State

*/

#include "EvaluatorIndependentFromFile.hh"
#include "EvaluatorFromFile_Helpers.hh"
#include "Function.hh"
#include "FunctionFactory.hh"

namespace Amanzi {

// ---------------------------------------------------------------------------
// Constructor
// ---------------------------------------------------------------------------
EvaluatorIndependentFromFile::EvaluatorIndependentFromFile(Teuchos::ParameterList& plist)
  : EvaluatorIndependent<CompositeVector, CompositeVectorSpace>(plist),
    filename_(plist.get<std::string>("filename")),
    meshname_(plist.get<std::string>("domain name", "domain")),
    compname_(plist.get<std::string>("component name", "cell")),
    varname_(plist.get<std::string>("variable name")),
    ndofs_(plist.get<int>("number of dofs", 1)),
    page_size_(plist.get<int>("page size", 2)),
    checkpoint_file_(plist.get<bool>("checkpoint file", false))
{
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

  // requirements on vector data
  auto& space = S.Require<CompositeVector, CompositeVectorSpace>(my_key_, my_tag_, my_key_);
  space.SetMesh(S.GetMesh(meshname_))->AddComponent(compname_, locname_, ndofs_);

  interpolator_ = Teuchos::rcp(new EvaluatorFromFile_Helpers::FileTimeInterpolator(
    filename_, varname_, compname_, ndofs_, checkpoint_file_, page_size_));
  interpolator_->setup(space, temporally_variable_);
}


// ---------------------------------------------------------------------------
// Update the value in the state.
// ---------------------------------------------------------------------------
void
EvaluatorIndependentFromFile::Update_(State& S)
{
  CompositeVector& cv = S.GetW<CompositeVector>(my_key_, my_tag_, my_key_);

  double t = S.get_time();
  if (time_func_ != Teuchos::null) {
    std::vector<double> point(1, t);
    t = (*time_func_)(point);
  }

  cv = interpolator_->interpolate(t);

  if (locname_ == AmanziMesh::Entity_kind::CELL &&
      (cv.HasComponent("boundary_face") || cv.HasComponent("face")))
    DeriveFaceValuesFromCellValues(cv);
}


} // namespace Amanzi
