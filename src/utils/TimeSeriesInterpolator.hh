/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

//! TimeSeriesInterpolator: pages a time series and interpolates in time.
/*
  Loading a range of time slices is done through an injected load_page
  callable; this class has zero knowledge of Input/FileHDF5/FileNetCDF, or
  of what the payload (every dimension of ViewType beyond dimension 0)
  actually represents.

  ViewType must be a contiguous, LayoutRight Kokkos::View of rank N+1:
  dimension 0 is the page's time-slot axis, every other dimension opaque
  payload.  page_storage and result_storage share this same type and
  payload shape (dimensions 1..rank-1); only dimension 0's extent differs
  -- page_size for page_storage, 1 for result_storage.  interpolate()
  returns a reference to result_storage itself, so callers index it as
  view(0, ...), not view(...).
*/
#pragma once

#include <algorithm>
#include <functional>
#include <vector>

#include "dbc.hh"

namespace Amanzi {

template <typename ViewType>
class TimeSeriesInterpolator {
 public:
  enum class Form_kind { LINEAR, CONSTANT };

  // load_page(start_time_index, count, dest_block) must fill page_storage's
  // slots [dest_block, dest_block + count) with the file's time indices
  // [start_time_index, start_time_index + count).
  using LoadPageFn = std::function<void(int, int, int)>;

  TimeSeriesInterpolator(std::vector<double> times,
                         int page_size,
                         LoadPageFn load_page,
                         ViewType page_storage,
                         ViewType result_storage,
                         Form_kind form = Form_kind::LINEAR)
    : times_(std::move(times)),
      page_size_(page_size),
      load_page_(std::move(load_page)),
      page_(std::move(page_storage)),
      result_(std::move(result_storage)),
      form_(form),
      page_start_(0),
      n_resident_(0)
  {
    AMANZI_ASSERT(page_size_ >= 2);
    AMANZI_ASSERT(!times_.empty());

    n_resident_ = std::min(page_size_, (int)times_.size());
    load_page_(0, n_resident_, 0);
  }

  // Ensures the resident page covers t, blends (or picks) the bracketing
  // slot(s) per form_, and returns the interpolated result.  Valid until
  // the next call to interpolate().
  const ViewType& interpolate(double t)
  {
    ensureWindow_(t);

    double t0 = times_[page_start_];
    if (n_resident_ == 1 || t <= t0) {
      // either at/before the first resident slice, or only one slice is
      // resident (past the end of the file, or a temporally-constant
      // series) -- hold constant.
      copySlot_(0, result_);
      return result_;
    }

    // find the resident slot pair straddling t: the largest slot such that
    // times_[page_start_ + slot] <= t.
    int slot = 0;
    while (slot + 1 < n_resident_ && t > times_[page_start_ + slot + 1]) ++slot;

    if (slot + 1 == n_resident_) {
      // t is beyond every resident time, and ensureWindow_ could not page
      // any further forward (the last resident slot is already the
      // series' last time index) -- hold the last resident slot constant
      // rather than reading times_[page_start_ + slot + 1], which may be
      // past the end of times_ itself.
      copySlot_(slot, result_);
      return result_;
    }

    // CONSTANT is a discontinuous form: a t that has drifted a hair short
    // of times_[slot+1] (e.g. a driver accumulating many small steps)
    // must still land on slot+1's value, not slot's -- 1.e-6 is chosen so
    // that dt in seconds is bigger than machine precision of a double on
    // 100 years in seconds.
    if (form_ == Form_kind::CONSTANT) {
      if (t + 1.e-6 > times_[page_start_ + slot + 1]) {
        copySlot_(slot + 1, result_);
      } else {
        copySlot_(slot, result_);
      }
      return result_;
    }

    if (t >= times_[page_start_ + slot + 1]) {
      // at the slot at slot+1 exactly -- return it rather than blending.
      copySlot_(slot + 1, result_);
      return result_;
    }

    double t_before = times_[page_start_ + slot];
    double t_after = times_[page_start_ + slot + 1];
    blendSlots_(t, t_before, t_after, slot, slot + 1, result_);
    return result_;
  }

 private:
  // Copies page slot `slot`'s payload into dest.
  void copySlot_(int slot, ViewType& dest)
  {
    const auto* src = page_.data() + slot * payloadSize_();
    auto* dst = dest.data();
    std::copy(src, src + payloadSize_(), dst);
  }

  // Linearly blends resident slots `slot_before` and `slot_after` at time
  // t, writing the result into dest.
  void blendSlots_(double t,
                   double t_before,
                   double t_after,
                   int slot_before,
                   int slot_after,
                   ViewType& dest)
  {
    AMANZI_ASSERT(t_after > t_before);
    double coef = (t - t_before) / (t_after - t_before);

    const auto* before = page_.data() + slot_before * payloadSize_();
    const auto* after = page_.data() + slot_after * payloadSize_();
    auto* dst = dest.data();
    std::size_t n = payloadSize_();
    for (std::size_t i = 0; i != n; ++i) dst[i] = (1 - coef) * before[i] + coef * after[i];
  }

  // Number of payload elements per page slot (i.e. page_.size() /
  // page_size_, the product of dimensions 1..rank-1).
  std::size_t payloadSize_() const { return page_.size() / page_size_; }

  // Ensures the resident page covers the interval containing t, rolling or
  // fully reloading the page as needed.
  void ensureWindow_(double t)
  {
    // degenerate case: at most one time slice exists -- nothing to page,
    // interpolate() always returns slot 0.
    if ((int)times_.size() <= 1) return;

    if (t < times_[page_start_] && page_start_ != 0) {
      // backward jump: restart by reloading the page from the beginning.
      page_start_ = 0;
      n_resident_ = std::min(page_size_, (int)times_.size());
      load_page_(0, n_resident_, 0);
      return;
    }

    // While t is within the outer bracket of the resident page, there is
    // nothing to do -- no I/O.  Only once t moves past the last resident
    // time do we need to page forward.
    while (t > times_[page_start_ + n_resident_ - 1] &&
           page_start_ + n_resident_ < (int)times_.size()) {
      if (n_resident_ < page_size_) {
        // page not yet full (near the start of a short series) -- just
        // load the next slice into the next free block.  No shifting
        // needed.
        load_page_(page_start_ + n_resident_, 1, n_resident_);
        n_resident_ += 1;
        continue;
      }

      // the page is full and t has moved past its far edge.  If t is
      // still within reach of paging forward by one page (i.e. within the
      // *next* page's span), do the cheap page-forward: keep the last
      // resident slice as the new first slice (it already brackets the
      // next interval), and load the remaining page_size_ - 1 new slices
      // in one range read.  Otherwise (a jump larger than a full page),
      // there is nothing reusable -- reload fresh starting near t.
      if (t > times_[std::min((int)times_.size() - 1, page_start_ + 2 * page_size_ - 2)]) {
        int new_start = page_start_;
        while (new_start + 1 < (int)times_.size() && times_[new_start + 1] <= t) ++new_start;
        page_start_ = new_start;
        n_resident_ = std::min(page_size_, (int)times_.size() - page_start_);
        load_page_(page_start_, n_resident_, 0);
        return;
      }

      // page forward: copy the last resident slice (slot n_resident_ - 1)
      // into slot 0 -- in-memory only, no file I/O -- then load the
      // remaining n_resident_ - 1 slices into slots [1, n_resident_) in
      // one range read.
      {
        auto* data = page_.data();
        std::size_t n = payloadSize_();
        std::copy(data + (n_resident_ - 1) * n, data + n_resident_ * n, data);
      }
      int new_page_start = page_start_ + n_resident_ - 1;
      int new_n_resident = std::min(page_size_, (int)times_.size() - new_page_start);
      if (new_n_resident > 1) load_page_(new_page_start + 1, new_n_resident - 1, 1);
      page_start_ = new_page_start;
      n_resident_ = new_n_resident;
    }
  }

  std::vector<double> times_;
  int page_size_;
  LoadPageFn load_page_;
  ViewType page_;
  ViewType result_;
  Form_kind form_;
  int page_start_; // index into times_ of the page's block 0
  int n_resident_; // number of time slices currently resident (<= page_size_)
};

} // namespace Amanzi
