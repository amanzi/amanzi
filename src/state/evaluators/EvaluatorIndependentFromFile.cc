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

#include "Reader.hh"

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
  S.Require<CompositeVector, CompositeVectorSpace>(my_key_, my_tag_, my_key_)
    .SetMesh(S.GetMesh(meshname_))
    ->AddComponent(compname_, locname_, ndofs_);

  // load times, ensure file is valid
  // if there exists no times, default value is set to +infinity
  auto reader = createReader(filename_);
  times_.clear();
  if (temporally_variable_) {
    try {
      Teuchos::Array<double> times;
      reader->read("/time", times);
      times_ = times.toVector();
    } catch (...) {
      std::stringstream messagestream;
      messagestream << "Variable " << my_key_ << " is defined as a field changing in time.\n"
                    << " Dataset /time is not provided in file " << filename_ << "\n";
      Errors::Message message(messagestream.str());
      Exceptions::amanzi_throw(message);
    }
  } else {
    times_.push_back(1e+99);
  }

  // check for increasing times
  for (int j = 1; j < times_.size(); ++j) {
    if (times_[j] <= times_[j - 1]) {
      Errors::Message m;
      m << "IndependentVariable from file: times values are not strictly "
           "increasing";
      throw(m);
    }
  }

  current_interval_ = -1;
  t_before_ = -1.0e+99;
  t_after_ = times_[0];
}


// ---------------------------------------------------------------------------
// Update the value in the state.
// ---------------------------------------------------------------------------
void
EvaluatorIndependentFromFile::Update_(State& S)
{
  CompositeVector& cv = S.GetW<CompositeVector>(my_key_, my_tag_, my_key_);

  if (!computed_once_) {
    val_after_ = Teuchos::rcp(new CompositeVector(cv));
    EvaluatorFromFile_Helpers::LoadFile(0, *this);
  }

  double t = S.get_time();
  if (time_func_ != Teuchos::null) {
    std::vector<double> point(1, t);
    t = (*time_func_)(point);
  }

  // check if we are before the current interval
  if (t < t_before_) {
    // this should only be the case if we are somehow composing this function
    // with a time function that is not monotonic, i.e. doing a cyclic steady
    // state to repeat a year, and we have gone past the cycle.  Restart the
    // interval.
    t_before_ = -1.0e+99;
    t_after_ = times_[0];
    current_interval_ = -1;
    EvaluatorFromFile_Helpers::LoadFile(0, *this);
  }

  // update by time interpolation
  EvaluatorFromFile_Helpers::UpdateTimeInterpolation(t, *this, cv);

  if (locname_ == AmanziMesh::Entity_kind::CELL &&
      (cv.HasComponent("boundary_face") || cv.HasComponent("face")))
    DeriveFaceValuesFromCellValues(cv);
}


} // namespace Amanzi
