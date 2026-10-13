/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

//! Tests TimeSeriesInterpolator's paging/windowing and LINEAR/CONSTANT
//! interpolation, against a toy (non-file-backed) series.

#include "UnitTest++.h"

#include <vector>

#include "Kokkos_Core.hpp"

#include "TimeSeriesInterpolator.hh"

using namespace Amanzi;

using ViewType = Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::HostSpace>;

// A toy series: time index k's payload is the 2-vector (k, k + 0.5), so
// LINEAR-interpolated values and the number of underlying slices loaded
// are both trivial to check by hand.  load_page_calls counts invocations
// of the fake load_page callable, to confirm paging only hits "disk" when
// it must.
struct ToySeries {
  int load_page_calls = 0;

  void load_page(ViewType page, int start, int count, int dest_block)
  {
    ++load_page_calls;
    for (int i = 0; i != count; ++i) {
      page(dest_block + i, 0) = start + i;
      page(dest_block + i, 1) = start + i + 0.5;
    }
  }
};


TEST(time_series_interpolator_linear_basic)
{
  // times = 0, 1, 2, 3, 4 ; page_size = 2
  std::vector<double> times = { 0., 1., 2., 3., 4. };
  int page_size = 2;
  ToySeries series;
  ViewType page("page", page_size, 2);
  ViewType result("result", 1, 2);

  TimeSeriesInterpolator<ViewType> interp(
    times,
    page_size,
    [&](int start, int count, int dest_block) { series.load_page(page, start, count, dest_block); },
    page,
    result);
  CHECK_EQUAL(1, series.load_page_calls);

  // t before the first time -- held at slice 0.
  auto v = interp.interpolate(-1.);
  CHECK_CLOSE(0., v(0, 0), 1.e-12);
  CHECK_CLOSE(0.5, v(0, 1), 1.e-12);

  // t exactly at a resident time -- returns that slice exactly, no blend.
  v = interp.interpolate(1.);
  CHECK_CLOSE(1., v(0, 0), 1.e-12);
  CHECK_CLOSE(1.5, v(0, 1), 1.e-12);

  // t strictly between two resident times -- linear blend.
  v = interp.interpolate(0.25);
  CHECK_CLOSE(0.25, v(0, 0), 1.e-12);
  CHECK_CLOSE(0.75, v(0, 1), 1.e-12);

  // Paging forward one slice at a time while the page isn't yet full.
  v = interp.interpolate(1.5);
  CHECK_CLOSE(1.5, v(0, 0), 1.e-12);
  CHECK_CLOSE(2.0, v(0, 1), 1.e-12);

  // t beyond the last time -- held at the last slice, no further I/O past
  // the end of the series.
  v = interp.interpolate(100.);
  CHECK_CLOSE(4., v(0, 0), 1.e-12);
  CHECK_CLOSE(4.5, v(0, 1), 1.e-12);
}


TEST(time_series_interpolator_linear_paging_forward)
{
  // A longer series (times 0..9), forcing the full page-forward (not
  // just-not-yet-full) path in ensureWindow_, and counting loads to
  // confirm I/O only happens when the resident window must move.
  // page_size = 3, so after construction the resident window is slots
  // {0,1,2} covering times [0,2].
  std::vector<double> times;
  for (int i = 0; i != 10; ++i) times.push_back((double)i);
  int page_size = 3;
  ToySeries series;
  ViewType page("page", page_size, 2);
  ViewType result("result", 1, 2);

  TimeSeriesInterpolator<ViewType> interp(
    times,
    page_size,
    [&](int start, int count, int dest_block) { series.load_page(page, start, count, dest_block); },
    page,
    result);
  CHECK_EQUAL(1, series.load_page_calls);

  // t within the resident window [0,2] -- no I/O.
  int calls_before = series.load_page_calls;
  auto v = interp.interpolate(1.5);
  CHECK_EQUAL(calls_before, series.load_page_calls);
  CHECK_CLOSE(1.5, v(0, 0), 1.e-12);

  // t = 2.5, past the window's far edge (2) but within reach of the cheap
  // page-forward (next page spans up to time 4) -- one load_page call,
  // window becomes {2,3,4} covering times [2,4].
  calls_before = series.load_page_calls;
  v = interp.interpolate(2.5);
  CHECK_EQUAL(calls_before + 1, series.load_page_calls);
  CHECK_CLOSE(2.5, v(0, 0), 1.e-12);

  // t = 3.5 is still within the new window [2,4] -- no I/O.
  calls_before = series.load_page_calls;
  v = interp.interpolate(3.5);
  CHECK_EQUAL(calls_before, series.load_page_calls);
  CHECK_CLOSE(3.5, v(0, 0), 1.e-12);

  // a big forward jump, beyond reach of the cheap page-forward, skips it
  // and reloads fresh starting at the bracketing time -- one load_page
  // call, window becomes {8,9} (times_.size() - 1 == 9 is the last index).
  calls_before = series.load_page_calls;
  v = interp.interpolate(8.5);
  CHECK_EQUAL(calls_before + 1, series.load_page_calls);
  CHECK_CLOSE(8.5, v(0, 0), 1.e-12);

  // a backward jump restarts from the beginning -- one load_page call.
  calls_before = series.load_page_calls;
  v = interp.interpolate(0.5);
  CHECK_EQUAL(calls_before + 1, series.load_page_calls);
  CHECK_CLOSE(0.5, v(0, 0), 1.e-12);
}


TEST(time_series_interpolator_constant)
{
  // times = 0, 1, 3 -- mirrors FunctionBilinearAndTime's test fixture
  // (functions/test/functions_all.cc's bilinear_and_time_test), where
  // CONSTANT holds the left (earlier) endpoint strictly inside an
  // interval, but an exact landing on a resident time (including via the
  // 1.e-6 roundoff nudge) takes that time's own value.
  std::vector<double> times = { 0., 1., 3. };
  int page_size = 3;
  ToySeries series;
  ViewType page("page", page_size, 2);
  ViewType result("result", 1, 2);

  TimeSeriesInterpolator<ViewType> interp(
    times,
    page_size,
    [&](int start, int count, int dest_block) { series.load_page(page, start, count, dest_block); },
    page,
    result,
    TimeSeriesInterpolator<ViewType>::Form_kind::CONSTANT);

  // exactly on time 1 -- takes slice 1's value, not slice 0's.
  auto v = interp.interpolate(1.0);
  CHECK_CLOSE(1., v(0, 0), 1.e-12);

  // strictly inside (1, 3) -- holds slice 1 (the left/earlier endpoint).
  v = interp.interpolate(2.0);
  CHECK_CLOSE(1., v(0, 0), 1.e-12);
  v = interp.interpolate(2.9);
  CHECK_CLOSE(1., v(0, 0), 1.e-12);

  // exactly on time 3 (the last resident time) -- takes slice 2's value.
  v = interp.interpolate(3.0);
  CHECK_CLOSE(2., v(0, 0), 1.e-12);

  // a hair short of time 3, within the 1.e-6 roundoff tolerance -- still
  // takes slice 2's value, not slice 1's.
  v = interp.interpolate(3.0 - 1.e-9);
  CHECK_CLOSE(2., v(0, 0), 1.e-12);
}
