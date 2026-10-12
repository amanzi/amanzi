/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

//! Tests TimeSeriesInterpolatorCV (EvaluatorFromFile_Helpers.hh).
/*
  Exercises TimeSeriesInterpolatorCV directly against an HDF5 file
  (bypassing the EvaluatorIndependentFromFile/State machinery, which is
  not needed -- TimeSeriesInterpolatorCV only needs a CompositeSpace), in
  particular the case where the time series length is not an even
  multiple of page_size -- i.e. the last resident page is short.  This
  exercises TimeSeriesInterpolator::ensureWindow_'s page-forward logic
  once page_start_ + n_resident_ reaches times_.size(): n_resident_ should
  shrink to the remainder rather than over-reading past the end of the
  file or under-filling silently.
*/

#include "UnitTest++.h"

#include "AmanziComm.hh"
#include "AmanziTypes.hh"
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "CompositeVectorSpace.hh"
#include "FileHDF5.hh"
#include "EvaluatorFromFile_Helpers.hh"

SUITE(EVALUATOR_FROM_FILE)
{
  using namespace Amanzi;

  // Writes an HDF5 file with ntimes time slices of a single-DoF, cell-based
  // field "myvar" on the given mesh's "cell" component, following
  // TimeSeriesInterpolatorCV/InputHDF5::readRange's on-disk convention:
  // /time (1D, length ntimes) and /myvar.cell.0/<time_index> (1D, one per
  // time index). value(t, i) gives the value at time index t, local cell i.
  template <typename ValueFn>
  void writeTimeSeries(const Comm_ptr_type& comm,
                       const std::string& filename,
                       const std::vector<double>& times,
                       int nlocal,
                       ValueFn value)
  {
    FileHDF5 out(comm, filename, FILE_CREATE);
    Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> time_view("time_view",
                                                                            times.size());
    for (std::size_t t = 0; t != times.size(); ++t) time_view(t) = times[t];
    out.writeLocalView("time", "/", time_view);

    for (int t = 0; t != (int)times.size(); ++t) {
      Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> vals("vals", nlocal);
      for (int i = 0; i != nlocal; ++i) vals(i) = value(t, i);
      out.writeDistributedView(std::to_string(t), "myvar.cell.0/", vals);
    }
  }


  // Builds a CompositeSpace for a single-DoF, cell-based "cell" component on
  // a small serial mesh -- everything TimeSeriesInterpolatorCV needs.
  Teuchos::RCP<const CompositeSpace> createCellSpace(
    const Teuchos::RCP<const AmanziMesh::Mesh>& mesh)
  {
    CompositeVectorSpace space;
    space.SetMesh(mesh);
    space.AddComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    return space.CreateSpace();
  }


  TEST(PAGE_SIZE_DOES_NOT_EVENLY_DIVIDE_TIMES)
  {
    auto comm = getDefaultComm();
    if (comm->getSize() > 1) return; // serial-only, mesh/comm details not under test

    AmanziMesh::MeshFactory fac(comm);
    fac.set_preference(AmanziMesh::Preference{ AmanziMesh::Framework::SIMPLE });
    auto mesh = fac.create(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2, 2, 2);
    int nlocal = mesh->getMap(AmanziMesh::Entity_kind::CELL, false)->getLocalNumElements();

    // 7 time slices, page_size 3 -- 7 % 3 != 0, so the last page is short
    // (page_start_=6, n_resident_=1 at the very end).
    std::vector<double> times{ 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0 };
    auto value = [&](int t, int i) { return 100 * t + i; };

    std::string filename = "test_state_evaluator_from_file_remainder.h5";
    writeTimeSeries(comm, filename, times, nlocal, value);

    auto space = createCellSpace(mesh);

    EvaluatorFromFile_Helpers::TimeSeriesInterpolatorCV interp(
      filename, "myvar", space, /*temporally_variable=*/true, /*page_size=*/3);
    CompositeVector cv(space);

    // walk every exact time slice, including the last one, and then well
    // past the end -- the latter must hold the final slice constant rather
    // than crashing or reading garbage.
    for (int t = 0; t != (int)times.size(); ++t) {
      interp.interpolate(times[t], cv);
      auto v = cv.viewComponent<MemSpace_kind::HOST>("cell", false);
      for (int i = 0; i != nlocal; ++i) CHECK_CLOSE(value(t, i), v(i, 0), 1.e-10);
    }

    {
      interp.interpolate(times.back() + 10.0, cv);
      auto v = cv.viewComponent<MemSpace_kind::HOST>("cell", false);
      for (int i = 0; i != nlocal; ++i) {
        CHECK_CLOSE(value((int)times.size() - 1, i), v(i, 0), 1.e-10);
      }
    }

    // and a value interpolated between the last two slices, to confirm the
    // short final page still blends correctly (not just holds constant).
    {
      interp.interpolate(5.5, cv);
      auto v = cv.viewComponent<MemSpace_kind::HOST>("cell", false);
      for (int i = 0; i != nlocal; ++i) {
        double expected = 0.5 * value(5, i) + 0.5 * value(6, i);
        CHECK_CLOSE(expected, v(i, 0), 1.e-10);
      }
    }
  }
}
