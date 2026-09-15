/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon, Bo Gao
*/

/*
  Share the following for CompositeVector and TensorVector-based evaluators

  Load data from HDF5 file
  Interpolate between two time steps
  Update the value in the state
  CopyVectorToTensorVector
  
*/


#pragma once

#include "TensorVector.hh"


namespace Amanzi {
namespace EvaluatorFromFile_Helpers {

template<typename EvaluatorType>
void LoadFile(int i, EvaluatorType& eval) 
{
  // allocate data
  if (eval.val_after_ == Teuchos::null) {
    AMANZI_ASSERT(eval.val_before_ != Teuchos::null);
    eval.val_after_ = Teuchos::rcp(new CompositeVector(*eval.val_before_));
  }

  // open the file
  Teuchos::RCP<Amanzi::HDF5_MPI> file_input =
    Teuchos::rcp(new Amanzi::HDF5_MPI(eval.val_after_->Comm(), eval.filename_));
  file_input->open_h5file();
  
  // load the data
  Epetra_MultiVector& vec = *eval.val_after_->ViewComponent(eval.compname_, false);
  for (int j = 0; j != eval.ndofs_; ++j) {
    std::stringstream varname;
    varname << eval.varname_ << "." << eval.compname_ << "." << j;
    if (!eval.checkpoint_file_) {
      varname << "//" << i;
    }
    file_input->readData(*vec(j), varname.str());
  }
  
  // close file
  file_input->close_h5file();
}


inline void Interpolate(double time, double t_before, double t_after,
                        const CompositeVector& val_before,
                        const CompositeVector& val_after,
                        CompositeVector& v)
{
  AMANZI_ASSERT(t_before >= 0.0);
  AMANZI_ASSERT(t_after >= 0.0);
  AMANZI_ASSERT(t_after >= time);
  AMANZI_ASSERT(time >= t_before);
  AMANZI_ASSERT(t_after > t_before);
 
  double coef = (time - t_before) / (t_after - t_before);
  v = val_before;
  v.Update(coef, val_after, 1 - coef);
}


template<typename EvaluatorType>
void UpdateTimeInterpolation(double t, EvaluatorType& eval, CompositeVector& cv)
{
  cv = *eval.val_after_;

  // Determine where we are relative to the currently stored interval
  if (t < eval.t_before_) {
    // should never be possible thanks to the previous check
    AMANZI_ASSERT(0);
  } else if (t == eval.t_before_) {
    // at the start of the interval
    AMANZI_ASSERT(eval.val_before_ != Teuchos::null);
    cv = *eval.val_before_;

  } else if (t < eval.t_after_) {
    if (eval.t_before_ == std::numeric_limits<double>::lowest()) {
      // to the left of the first point
      AMANZI_ASSERT(eval.val_after_ != Teuchos::null);
      cv = *eval.val_after_;
    } else if (eval.val_after_ == Teuchos::null) {
      // to the right of the last point
      AMANZI_ASSERT(eval.val_before_ != Teuchos::null);
      cv = *eval.val_before_;
    } else {
      // in the interval, interpolate
      Interpolate(t, eval.t_before_, eval.t_after_, 
                  *eval.val_before_, *eval.val_after_, cv);
    }
  } else if (t == eval.t_after_) {
    // at the end of the interval
    AMANZI_ASSERT(eval.val_after_ != Teuchos::null);
    cv = *eval.val_after_;

  } else {
    // to the right of the interval -- advance the interval
    while (t > eval.t_after_) {
      eval.current_interval_++;
      eval.t_before_ = eval.t_after_;
      if (eval.current_interval_ + 1 == eval.times_.size()) {
        // at the end of data
        eval.t_after_ = 1.0e99;
        eval.val_before_ = eval.val_after_;
        eval.val_after_ = Teuchos::null;

        // copy the value
        cv = *eval.val_before_;

      } else {
        eval.t_after_ = eval.times_[eval.current_interval_ + 1];

        // swap the pointers
        std::swap(eval.val_before_, eval.val_after_);

        // load the new data
        LoadFile(eval.current_interval_ + 1, eval);
 
        // now we are in the interval, interpolate
        if (t == eval.t_after_) {
          cv = *eval.val_after_;
        } else if (t < eval.t_after_) {
          Interpolate(t, eval.t_before_, eval.t_after_,
                      *eval.val_before_, *eval.val_after_, cv);
        }
      }
    }
  }
}


inline void CopyVectorToTensorVector(const Epetra_MultiVector& v, int j, TensorVector& tv)
{
  AMANZI_ASSERT(v.MyLength() == tv.size());

  unsigned int ni = v.MyLength();
  unsigned int ndofs = v.NumVectors();
  unsigned int space_dim = tv.dim;

  if (ndofs == 1) { // isotropic
    for (unsigned int i = 0; i != ni; ++i) tv[i + j](0, 0) = v[0][i];

  } else if (ndofs == 2 && space_dim == 3) {
    // horizontal and vertical perms
    for (int i = 0; i != ni; ++i) {
      tv[i + j](0, 0) = v[0][i];
      tv[i + j](1, 1) = v[0][i];
      tv[i + j](2, 2) = v[1][i];
    }

  } else if (ndofs >= space_dim) {
    // diagonal tensor
    for (unsigned int dim = 0; dim != space_dim; ++dim) {
      for (unsigned int i = 0; i != ni; ++i) {
        tv[i + j](dim, dim) = v[dim][i];
      }
    }

    if (ndofs > space_dim) {
      // full tensor
      if (ndofs == 3) { // 2D
        for (unsigned int i = 0; i != ni; ++i) {
          tv[i + j](0, 1) = tv[i + j](1, 0) = v[2][i];
        }

      } else if (ndofs == 6) { // 3D
        for (unsigned int i = 0; i != ni; ++i) {
          tv[i + j](0, 1) = tv[i + j](1, 0) = v[3][i]; // xy & yx
          tv[i + j](0, 2) = tv[i + j](2, 0) = v[4][i]; // xz & zx
          tv[i + j](1, 2) = tv[i + j](2, 1) = v[5][i]; // yz & zy
        }
      } else if (ndofs == 9) { // 3D full tensor
        for (unsigned int i = 0; i != ni; ++i) {
          tv[i + j](0, 0) = v[0][i];
          tv[i + j](0, 1) = v[1][i];
          tv[i + j](0, 2) = v[2][i];
          tv[i + j](1, 0) = v[3][i];
          tv[i + j](1, 1) = v[4][i];
          tv[i + j](1, 2) = v[5][i];
          tv[i + j](2, 0) = v[6][i];
          tv[i + j](2, 1) = v[7][i];
          tv[i + j](2, 2) = v[8][i];
        }
      } else {
        AMANZI_ASSERT(0);
      }
    }

  } else {
    // ERROR -- unknown perm type
    AMANZI_ASSERT(0);
  }
  
}


} // namespace EvaluatorFromFile_Helpers
} // namespace Amanzi
