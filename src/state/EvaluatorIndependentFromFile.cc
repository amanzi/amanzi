/*
  State

  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon
*/

#include "HDF5Reader.hh"

#include "EvaluatorIndependentFromFile.hh"
#include "Function.hh"
#include "FunctionFactory.hh"

namespace Amanzi {

// ---------------------------------------------------------------------------
// Constructor
// ---------------------------------------------------------------------------
EvaluatorIndependentFromFile::EvaluatorIndependentFromFile(
    Teuchos::ParameterList& plist)
    : EvaluatorIndependent<CompositeVector, CompositeVectorSpace>(plist),
      filename_(plist.get<std::string>("filename")),
      meshname_(plist.get<std::string>("domain name", "domain")),
      varname_(plist.get<std::string>("variable name")),
      compname_(plist.get<std::string>("component name", "cell")),
      locname_(plist.get<std::string>("mesh entity", "cell")),
      ndofs_(plist.get<int>("number of dofs", 1)) {
  if (plist.isSublist("time function")) {
    FunctionFactory fac;
    time_func_ = Teuchos::rcp(fac.Create(plist.sublist("time function")));
  }
}


// ---------------------------------------------------------------------------
// Virtual Copy constructor
// ---------------------------------------------------------------------------
Teuchos::RCP<Evaluator> EvaluatorIndependentFromFile::Clone() const {
  return Teuchos::rcp(new EvaluatorIndependentFromFile(*this));
}


// ---------------------------------------------------------------------------
// Operator=
// ---------------------------------------------------------------------------
Evaluator& EvaluatorIndependentFromFile::operator=(const Evaluator& other) {
  if (this !=& other) {
    const EvaluatorIndependentFromFile *other_p =
        dynamic_cast<const EvaluatorIndependentFromFile *>(&other);
    AMANZI_ASSERT(other_p != NULL);
    *this = *other_p;
  }
  return *this;
}


EvaluatorIndependentFromFile& EvaluatorIndependentFromFile::
operator=(const EvaluatorIndependentFromFile& other) {
  if (this !=& other) {
    AMANZI_ASSERT(my_key_ == other.my_key_);
    requests_ = other.requests_;
  }
  return *this;
}


// ---------------------------------------------------------------------------
// Ensures that the function can provide for the vector's requirements.
// ---------------------------------------------------------------------------
void EvaluatorIndependentFromFile::EnsureCompatibility(State& S) {
  EvaluatorIndependent::EnsureCompatibility(S);

  // requirements on vector data
  if (locname_ == "cell") {
    S.Require<CompositeVector, CompositeVectorSpace>(my_key_, my_tag_, my_key_)
        .SetMesh(S.GetMesh(meshname_))
        ->AddComponent(compname_, AmanziMesh::CELL, ndofs_);
  } else if (locname_ == "face") {
    S.Require<CompositeVector, CompositeVectorSpace>(my_key_, my_tag_, my_key_)
        .SetMesh(S.GetMesh(meshname_))
        ->AddComponent(compname_, AmanziMesh::FACE, ndofs_);
  } else if (locname_ == "boundary_face") {
    S.Require<CompositeVector, CompositeVectorSpace>(my_key_, my_tag_, my_key_)
        .SetMesh(S.GetMesh(meshname_))
        ->AddComponent(compname_, AmanziMesh::BOUNDARY_FACE, ndofs_);
  } else {
    Errors::Message m;
    m << "IndependentVariableFromFile: invalid location name: \"" << locname_
      << "\"";
    throw(m);
  }

  // load times, ensure file is valid
  // if there exists no times, default value is set to +infinity
  HDF5Reader reader(filename_);
  times_.clear();
  if (temporally_variable_) {
    try {
      reader.ReadData("/time", times_);
    } catch (...) {
      std::stringstream messagestream;
      messagestream << "Variable "<< my_key_ << " is defined as a field changing in time.\n"
                    << " Dataset /time is not provided in file " << filename_<<"\n";
      Errors::Message message(messagestream.str());
      Exceptions::amanzi_throw(message);
    }
  } else{
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
void EvaluatorIndependentFromFile::Update_(State& S)
{
  CompositeVector& cv = S.GetW<CompositeVector>(my_key_, my_tag_, my_key_);

  if (!computed_once_) {
    val_after_ = Teuchos::rcp(new CompositeVector(cv));
    LoadFile_(0);
  }

  double t = S.time();
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
    LoadFile_(0);
  }

  // determine where we are relative to the currently stored interval
  if (t < t_before_) {
    // should never be possible thanks to the previous check
    AMANZI_ASSERT(0);
  } else if (t == t_before_) {
    // at the start of the interval
    AMANZI_ASSERT(val_before_ != Teuchos::null);
    cv = *val_before_;

  } else if (t < t_after_) {
    if (t_before_ == -1.0e+99) {
      // to the left of the first point
      AMANZI_ASSERT(val_after_ != Teuchos::null);
      cv = *val_after_;
    } else if (val_after_ == Teuchos::null) {
      // to the right of the last point
      AMANZI_ASSERT(val_before_ != Teuchos::null);
      cv = *val_before_;
    } else {
      // in the interval, interpolate
      Interpolate_(t, cv);
    }
  } else if (t == t_after_) {
    // at the end of the interval
    AMANZI_ASSERT(val_after_ != Teuchos::null);
    cv = *val_after_;

  } else {
    // to the right of the interval -- advance the interval
    while (t > t_after_) {
      current_interval_++;
      t_before_ = t_after_;
      if (current_interval_ + 1 == times_.size()) {
        // at the end of data
        t_after_ = 1.0e99;
        val_before_ = val_after_;
        val_after_ = Teuchos::null;

        // copy the value
        cv = *val_before_;

      } else {
        t_after_ = times_[current_interval_ + 1];

        // swap the pointers
        std::swap(val_before_, val_after_);

        // load the new data
        LoadFile_(current_interval_ + 1);

        // now we are in the interval, interpolate
        if (t == t_after_) {
          cv = *val_after_;
        } else if (t < t_after_) {
          Interpolate_(t, cv);
        }
      }
    }
  }

  if (locname_ == "cell" &&
      (cv.HasComponent("boundary_face") || cv.HasComponent("face")))
    DeriveFaceValuesFromCellValues(cv);
}


// ---------------------------------------------------------------------------
//
// ---------------------------------------------------------------------------
void EvaluatorIndependentFromFile::LoadFile_(int i)
{
  // allocate data
  if (val_after_ == Teuchos::null) {
    AMANZI_ASSERT(val_before_ != Teuchos::null);
    val_after_ = Teuchos::rcp(new CompositeVector(*val_before_));
  }

  // open the file
  Teuchos::RCP<Amanzi::HDF5_MPI> file_input =
      Teuchos::rcp(new Amanzi::HDF5_MPI(val_after_->Comm(), filename_));
  file_input->open_h5file();

  // load the data
  Epetra_MultiVector& vec = *val_after_->ViewComponent(compname_, false);
  for (int j = 0; j != ndofs_; ++j) {
    std::stringstream varname;
    varname << varname_ << "." << locname_ << "." << j << "//" << i;
    file_input->readData(*vec(j), varname.str());
  }

  // close file
  file_input->close_h5file();
}


// ---------------------------------------------------------------------------
//
// ---------------------------------------------------------------------------
void EvaluatorIndependentFromFile::Interpolate_(double time, CompositeVector& v)
{
  AMANZI_ASSERT(t_before_ >= 0.0);
  AMANZI_ASSERT(t_after_ >= 0.0);
  AMANZI_ASSERT(t_after_ >= time);
  AMANZI_ASSERT(time >= t_before_);
  AMANZI_ASSERT(t_after_ > t_before_);
  AMANZI_ASSERT(val_before_ != Teuchos::null);
  AMANZI_ASSERT(val_after_ != Teuchos::null);

  double coef = (time - t_before_) / (t_after_ - t_before_);
  v = *val_before_;
  v.Update(coef, *val_after_, 1 - coef);
}

} // namespace Amanzi
