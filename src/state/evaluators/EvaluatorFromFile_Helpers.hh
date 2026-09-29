/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

//! Shared machinery for "from file" independent evaluators.
/*
  State

  FileTimeInterpolator owns reading a series of time slices of a
  CompositeVector-shaped variable from an HDF5 file and linearly
  interpolating between them in time.  It is shared by
  EvaluatorIndependentFromFile and EvaluatorIndependentTensorFromFile (the
  latter converts the resulting, per-timestep CompositeVector into a
  TensorVector after interpolate() returns).

  Within the file, data is expected to meet the following (HDF5) layout::

     /time : a 1D array of length NTIMES, providing the time in seconds.
     /variable_name.ENTITY.DOF  (group)

        /0 : a 1D array of length NENTITIES, providing the values for each
             entity at time /time[0]
        /1 : ...
        /NTIMES-1 : 1D array at time /time[NTIMES-1]

  or, if checkpoint_file is true, /variable_name.ENTITY.DOF is itself a
  single 1D array (no time index/group), and the variable is assumed
  constant in time.
*/

#pragma once

#include <limits>
#include <string>
#include <utility>
#include <vector>

#include "Teuchos_RCP.hpp"

#include "errors.hh"
#include "exceptions.hh"

#include "HDF5Reader.hh"
#include "FileHDF5.hh"
#include "CompositeVector.hh"
#include "CompositeVectorSpace.hh"

namespace Amanzi {
namespace EvaluatorFromFile_Helpers {

//
// Owns the state needed to read a time-series variable from an HDF5 file and
// linearly interpolate it in time, keeping a resident "page" of up to
// page_size_ consecutive time slices in memory to avoid reloading from file
// on every interval advance.
//
// The page is stored as a single CompositeVector whose components have
// ndofs_ * page_size_ columns (page_size_ contiguous blocks of ndofs_
// columns each), rather than as page_size_ separate CompositeVector
// objects.  This keeps the door open for reading a whole page in one
// bulk/hyperslab read in the future -- the only place that talks to the
// file is loadRange_(), which already takes a range of time indices, even
// though today it is implemented as a loop of single-time-index reads.
//
class FileTimeInterpolator {
 public:
  FileTimeInterpolator(std::string filename,
                       std::string varname,
                       std::string compname,
                       std::string locname,
                       int ndofs,
                       bool checkpoint_file,
                       int page_size = 2)
    : filename_(std::move(filename)),
      varname_(std::move(varname)),
      compname_(std::move(compname)),
      locname_(std::move(locname)),
      ndofs_(ndofs),
      checkpoint_file_(checkpoint_file),
      page_size_(page_size),
      page_start_(0),
      n_resident_(0)
  { AMANZI_ASSERT(page_size_ >= 2); }

  // Reads /time (if temporally_variable), validates strictly increasing
  // times, and loads the initial page starting at time index 0.  Must be
  // called once, after the CompositeVectorSpace defining the data's layout
  // is finalized.
  void setup(const Teuchos::RCP<const CompositeSpace>& space, bool temporally_variable)
  {
    space_ = space;
    temporally_variable_ = temporally_variable;

    if (temporally_variable_) {
      // NOTE: all ranks read this.  Could do on rank 0 and broadcast.
      HDF5Reader reader(filename_);
      reader.ReadData("/time", times_);

      for (int j = 1; j < times_.size(); ++j) {
        if (times_[j] <= times_[j - 1]) {
          Errors::Message m;
          m << "FileTimeInterpolator: times values in \"" << filename_
            << "\" are not strictly increasing";
          Exceptions::amanzi_throw(m);
        }
      }
    } else {
      times_.clear();
      times_.push_back(std::numeric_limits<double>::max());
    }

    // build the page's space: same components/locations/mesh as space_, but
    // ndofs_ * page_size_ columns per component instead of ndofs_.
    CompositeVectorSpace page_space;
    page_space.SetMesh(space_->getMesh());
    for (const auto& name : *space_) {
      page_space.AddComponent(name, space_->getLocation(name), ndofs_ * page_size_);
    }
    page_ = Teuchos::rcp(new CompositeVector(page_space.CreateSpace()));
    interpolated_ = Teuchos::rcp(new CompositeVector(space_));

    page_start_ = 0;
    n_resident_ = std::min(page_size_, (int)times_.size());
    loadRange_(0, n_resident_, 0);
  }

  // Returns the (possibly repaged) interpolated value at time t.  The
  // returned reference is valid until the next call to interpolate().
  const CompositeVector& interpolate(double t)
  {
    ensureWindow_(t);

    double t0 = times_[page_start_];
    if (n_resident_ == 1 || t <= t0) {
      // either at/before the first resident slice, or only one slice is
      // resident (past the end of the file, or the checkpoint_file_ /
      // temporally-constant case) -- hold constant.
      copyColumn_(0, *interpolated_);
      return *interpolated_;
    }

    // find the resident slot pair straddling t: the largest slot such that
    // times_[page_start_ + slot] <= t.
    int slot = 0;
    while (slot + 1 < n_resident_ && t > times_[page_start_ + slot + 1]) ++slot;

    if (t >= times_[page_start_ + slot + 1]) {
      // at (or, if this is the last resident slot pair and there is no more
      // data, beyond) the slot at slot+1 -- return it exactly rather than
      // blending.
      copyColumn_(slot + 1, *interpolated_);
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
  // though it is currently implemented as a loop of single-time-index
  // reads, so that a future bulk/hyperslab read can be plugged in here
  // without touching the rest of this class.
  void loadRange_(int start_time_index, int count, int dest_block)
  {
    FileHDF5 file_input(page_->getComm(), filename_, FILE_READONLY);

    MultiVector_type& vec = *page_->getComponent(compname_, false);
    // historically these files are written by component: dataset path is
    // "varname_/.locname_.dof" (a group per DoF), holding one dataset per
    // time index -- unless this is a checkpoint-style file, in which case
    // the group itself is the (single, constant-in-time) dataset.
    for (int c = 0; c != count; ++c) {
      int time_index = start_time_index + c;
      int block = dest_block + c;
      for (int j = 0; j != ndofs_; ++j) {
        std::string name = varname_ + "/" + "." + locname_ + "." + std::to_string(j);
        if (!checkpoint_file_) name += "/" + std::to_string(time_index);
        file_input.readVector(name, *vec.getVectorNonConst(block * ndofs_ + j));
      }
    }
  }

  // Copies page column-block `slot` into dest (dest has ndofs_ columns per
  // component, page_ has ndofs_ * page_size_).
  void copyColumn_(int slot, CompositeVector& dest)
  {
    for (const auto& name : dest) {
      MultiVector_type& dvec = *dest.getComponent(name, false);
      MultiVector_type& pvec = *page_->getComponent(name, false);
      for (int j = 0; j != ndofs_; ++j) {
        dvec.getVectorNonConst(j)->assign(*pvec.getVectorNonConst(slot * ndofs_ + j));
      }
    }
  }

  // Linearly blends resident slots `slot_before` and `slot_after` at time
  // t, writing the result into dest.
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
      MultiVector_type& dvec = *dest.getComponent(name, false);
      MultiVector_type& pvec = *page_->getComponent(name, false);
      for (int j = 0; j != ndofs_; ++j) {
        dvec.getVectorNonConst(j)->update(1 - coef,
                                          *pvec.getVectorNonConst(slot_before * ndofs_ + j),
                                          coef,
                                          *pvec.getVectorNonConst(slot_after * ndofs_ + j),
                                          0.0);
      }
    }
  }

  // Ensures the resident page covers the interval containing t, rolling or
  // fully reloading the page as needed.
  void ensureWindow_(double t)
  {
    // degenerate case: at most one time slice exists in the file (the
    // checkpoint_file_ / temporally-constant case, or a file with a single
    // time slice) -- nothing to page, interpolate() always returns slot 0.
    if ((int)times_.size() <= 1) return;

    if (t < times_[page_start_] && page_start_ != 0) {
      // backward jump: restart by reloading the page from the beginning.
      // E.g. a time function composed with a non-monotonic (e.g. cyclic)
      // time is going backward past the start of the resident window.  Note
      // the page_start_ != 0 guard -- if the page is already at the start,
      // there is nothing to reload; t before times_[0] is handled by
      // interpolate() holding the first slice constant, and without this
      // guard every such call would trigger a pointless reload.
      page_start_ = 0;
      n_resident_ = std::min(page_size_, (int)times_.size());
      loadRange_(0, n_resident_, 0);
      return;
    }

    // While t is within the outer bracket of the resident page, there is
    // nothing to do -- no I/O.  Only once t moves past the last resident
    // time do we need to page forward.
    while (t > times_[page_start_ + n_resident_ - 1] &&
           page_start_ + n_resident_ < (int)times_.size()) {
      if (n_resident_ < page_size_) {
        // page not yet full (near the start of a short file) -- just load
        // the next slice into the next free block.  No shifting needed.
        loadRange_(page_start_ + n_resident_, 1, n_resident_);
        n_resident_ += 1;
        continue;
      }

      // the page is full and t has moved past its far edge.  If t is still
      // within reach of paging forward by one page (i.e. within the *next*
      // page's span), do the cheap page-forward: keep the last resident
      // slice as the new first slice (it already brackets the next
      // interval), and load the remaining page_size_ - 1 new slices in one
      // range read.  Otherwise (a jump larger than a full page), there is
      // nothing reusable -- reload fresh starting near t.
      if (t > times_[std::min((int)times_.size() - 1, page_start_ + 2 * page_size_ - 2)]) {
        int new_start = page_start_;
        while (new_start + 1 < (int)times_.size() && times_[new_start + 1] <= t) ++new_start;
        page_start_ = new_start;
        n_resident_ = std::min(page_size_, (int)times_.size() - page_start_);
        loadRange_(page_start_, n_resident_, 0);
        return;
      }

      // page forward: copy the last resident slice (block n_resident_ - 1)
      // into block 0 -- pointer data only, no file I/O -- then load the
      // remaining n_resident_ - 1 slices into blocks [1, n_resident_) in one
      // range read.
      for (const auto& name : *page_) {
        MultiVector_type& pvec = *page_->getComponent(name, false);
        for (int j = 0; j != ndofs_; ++j) {
          pvec.getVectorNonConst(j)->assign(
            *pvec.getVectorNonConst((n_resident_ - 1) * ndofs_ + j));
        }
      }
      int new_page_start = page_start_ + n_resident_ - 1;
      int new_n_resident = std::min(page_size_, (int)times_.size() - new_page_start);
      if (new_n_resident > 1) loadRange_(new_page_start + 1, new_n_resident - 1, 1);
      page_start_ = new_page_start;
      n_resident_ = new_n_resident;
    }
  }

  std::string filename_, varname_, compname_, locname_;
  int ndofs_;
  bool checkpoint_file_;
  bool temporally_variable_;
  int page_size_;

  Teuchos::RCP<const CompositeSpace> space_;
  std::vector<double> times_;
  int page_start_; // index into times_ of the page's block 0
  int n_resident_; // number of time slices currently resident (<= page_size_)

  Teuchos::RCP<CompositeVector> page_;
  Teuchos::RCP<CompositeVector> interpolated_; // scratch, returned by interpolate()
};

} // namespace EvaluatorFromFile_Helpers
} // namespace Amanzi
