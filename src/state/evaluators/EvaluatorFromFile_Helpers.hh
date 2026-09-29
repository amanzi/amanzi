/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon, Bo Gao
*/

/*
  Shared machinery for CompositeVector and TensorVector-based "from file"
  evaluators:

  - FileTimeInterpolator: owns reading a series of time slices of a variable
    from an HDF5 file and interpolating between them in time.  Keeps a
    resident "page" of several consecutive time slices in memory to avoid
    reloading from file on every interval advance.
  - CopyVectorToTensorVector: stateless conversion of a CompositeVector's
    per-dof layout into a TensorVector.
*/


#pragma once

#include <string>
#include <vector>

#include "Teuchos_RCP.hpp"

#include "errors.hh"
#include "exceptions.hh"
#include "HDF5_MPI.hh"
#include "Reader.hh"
#include "CompositeVector.hh"
#include "CompositeVectorSpace.hh"
#include "TensorVector.hh"


namespace Amanzi {
namespace EvaluatorFromFile_Helpers {

//
// Owns the state needed to read a time-series variable from an HDF5 file and
// linearly interpolate it in time, keeping a resident "page" of up to
// page_size consecutive time slices in memory.
//
// The page is stored as a single CompositeVector whose components have
// ndofs * page_size columns (page_size contiguous blocks of ndofs columns
// each), rather than as page_size separate CompositeVector objects.  This
// keeps the door open for reading a whole page in one bulk/hyperslab read in
// the future -- the only place that talks to the file is loadRange_(), which
// already takes a range of time indices, even though today it is
// implemented as a loop of single-time-index reads.
//
class FileTimeInterpolator {
 public:
  FileTimeInterpolator(std::string filename,
                       std::string varname,
                       std::string compname,
                       int ndofs,
                       bool checkpoint_file,
                       int page_size = 2)
    : filename_(std::move(filename)),
      varname_(std::move(varname)),
      compname_(std::move(compname)),
      ndofs_(ndofs),
      checkpoint_file_(checkpoint_file),
      page_size_(page_size),
      page_start_(0)
  {
    AMANZI_ASSERT(page_size_ >= 2);
  }

  // Reads /time (if temporally_variable), validates strictly increasing
  // times, and loads the initial page starting at time index 0.
  void setup(const CompositeVectorSpace& space, bool temporally_variable)
  {
    auto reader = createReader(filename_);
    times_.clear();
    if (temporally_variable) {
      try {
        Teuchos::Array<double> times;
        reader->read("/time", times);
        times_ = times.toVector();
      } catch (...) {
        std::stringstream messagestream;
        messagestream << "FileTimeInterpolator: variable " << varname_
                      << " is defined as a field changing in time.\n"
                      << " Dataset /time is not provided in file " << filename_ << "\n";
        Errors::Message message(messagestream.str());
        Exceptions::amanzi_throw(message);
      }
    } else {
      times_.push_back(std::numeric_limits<double>::max());
    }

    for (int j = 1; j < times_.size(); ++j) {
      if (times_[j] <= times_[j - 1]) {
        Errors::Message m;
        m << "FileTimeInterpolator: times values are not strictly increasing";
        Exceptions::amanzi_throw(m);
      }
    }

    // build the page's space: same components/locations, but ndofs *
    // page_size columns per component instead of ndofs.
    CompositeVectorSpace page_space;
    page_space.SetMesh(space.Mesh());
    for (const auto& name : space) {
      page_space.AddComponent(name, space.Location(name), ndofs_ * page_size_);
    }
    page_ = Teuchos::rcp(new CompositeVector(page_space));
    interpolated_ = Teuchos::rcp(new CompositeVector(space));

    page_start_ = 0;
    n_resident_ = std::min(page_size_, (int)times_.size());
    loadRange_(0, n_resident_, 0);
  }

  // Returns the (possibly repaged) interpolated value at time t.
  const CompositeVector& interpolate(double t)
  {
    ensureWindow_(t);

    double t0 = times_[page_start_];
    if (n_resident_ == 1 || t <= t0) {
      // either at/before the first resident slice, or only one slice is
      // resident (past the end of the file) -- hold constant.
      copyColumn_(0, *interpolated_);
      return *interpolated_;
    }

    // find the resident slot straddling t
    int slot = 0;
    while (slot + 1 < n_resident_ && t > times_[page_start_ + slot + 1]) ++slot;

    if (slot + 1 >= n_resident_ || t >= times_[page_start_ + slot + 1]) {
      // at or beyond the last resident slice
      copyColumn_(n_resident_ - 1, *interpolated_);
      return *interpolated_;
    }

    double t_before = times_[page_start_ + slot];
    double t_after = times_[page_start_ + slot + 1];
    blendColumns_(t, t_before, t_after, slot, slot + 1, *interpolated_);
    return *interpolated_;
  }

 private:
  // Reads time indices [start_time_index, start_time_index + count) into
  // page column-blocks [dest_block, dest_block + count).  This is the one
  // place that talks to the file; it takes a *range* of time indices even
  // though it is currently implemented as a loop of single-time-index reads,
  // so that a future bulk/hyperslab read can be plugged in here without
  // touching the rest of this class.
  void loadRange_(int start_time_index, int count, int dest_block)
  {
    Teuchos::RCP<Amanzi::HDF5_MPI> file_input =
      Teuchos::rcp(new Amanzi::HDF5_MPI(page_->Comm(), filename_));
    file_input->open_h5file();

    Epetra_MultiVector& vec = *page_->ViewComponent(compname_, false);
    for (int c = 0; c != count; ++c) {
      int time_index = start_time_index + c;
      int block = dest_block + c;
      for (int j = 0; j != ndofs_; ++j) {
        std::stringstream varname;
        varname << varname_ << "." << compname_ << "." << j;
        if (!checkpoint_file_) {
          varname << "//" << time_index;
        }
        file_input->readData(*vec(block * ndofs_ + j), varname.str());
      }
    }

    file_input->close_h5file();
  }

  // Copies page column-block `slot` into dest (dest has ndofs_ columns per
  // component, page_ has ndofs_ * page_size_).
  void copyColumn_(int slot, CompositeVector& dest)
  {
    for (const auto& name : dest) {
      Epetra_MultiVector& dvec = *dest.ViewComponent(name, false);
      Epetra_MultiVector& pvec = *page_->ViewComponent(name, false);
      for (int j = 0; j != ndofs_; ++j) *dvec(j) = *pvec(slot * ndofs_ + j);
    }
  }

  // Linearly blends resident slots `slot_before` and `slot_after` at time t,
  // writing the result into dest.
  void blendColumns_(double t,
                     double t_before,
                     double t_after,
                     int slot_before,
                     int slot_after,
                     CompositeVector& dest)
  {
    AMANZI_ASSERT(t_after > t_before);
    double coef = (t - t_before) / (t_after - t_before);

    for (const auto& name : dest) {
      Epetra_MultiVector& dvec = *dest.ViewComponent(name, false);
      Epetra_MultiVector& pvec = *page_->ViewComponent(name, false);
      for (int j = 0; j != ndofs_; ++j) {
        dvec(j)->Update(
          1 - coef, *pvec(slot_before * ndofs_ + j), coef, *pvec(slot_after * ndofs_ + j), 0.0);
      }
    }
  }

  // Ensures the resident page covers the interval containing t, rolling or
  // fully reloading the page as needed.
  void ensureWindow_(double t)
  {
    if (t < times_[page_start_]) {
      // restart: reload the page from the beginning
      page_start_ = 0;
      n_resident_ = std::min(page_size_, (int)times_.size());
      loadRange_(0, n_resident_, 0);
      return;
    }

    // advance while t is beyond the last *resident* time and there is more
    // data to page in
    while (t > times_[page_start_ + n_resident_ - 1] &&
           page_start_ + n_resident_ < (int)times_.size()) {
      // how far are we jumping? if the jump spans past the whole current
      // page, there is nothing reusable -- reload fresh starting near t.
      // (single-step advance -- one interval at a time -- is the common
      // case and is handled below without any re-reading of resident data.)
      if (t > times_[std::min((int)times_.size() - 1, page_start_ + page_size_)]) {
        int new_start = page_start_;
        while (new_start + 1 < (int)times_.size() && times_[new_start + 1] <= t) ++new_start;
        page_start_ = new_start;
        n_resident_ = std::min(page_size_, (int)times_.size() - page_start_);
        loadRange_(page_start_, n_resident_, 0);
        return;
      }

      // single-step advance: reuse resident columns, load only the new one
      if (n_resident_ == page_size_) {
        // shift blocks [1, page_size_) down to [0, page_size_-1) -- pointer
        // data only, no file I/O.
        for (const auto& name : *page_) {
          Epetra_MultiVector& pvec = *page_->ViewComponent(name, false);
          for (int block = 1; block != page_size_; ++block) {
            for (int j = 0; j != ndofs_; ++j) {
              *pvec((block - 1) * ndofs_ + j) = *pvec(block * ndofs_ + j);
            }
          }
        }
        page_start_ += 1;
        loadRange_(page_start_ + page_size_ - 1, 1, page_size_ - 1);
      } else {
        // page not yet full (near the start of a short file) -- just load
        // the next slice into the next free block.
        loadRange_(page_start_ + n_resident_, 1, n_resident_);
        n_resident_ += 1;
      }
    }
  }

  std::string filename_, varname_, compname_;
  int ndofs_;
  bool checkpoint_file_;
  int page_size_;

  std::vector<double> times_;
  int page_start_; // index into times_ of the page's block 0
  int n_resident_; // number of time slices currently resident (<= page_size_)

  Teuchos::RCP<CompositeVector> page_;
  Teuchos::RCP<CompositeVector> interpolated_; // scratch, returned by interpolate()
};


inline void
CopyVectorToTensorVector(const Epetra_MultiVector& v, int j, TensorVector& tv)
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
