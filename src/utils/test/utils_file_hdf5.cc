/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

//! Tests FileHDF5's read/write API: Attribute, LocalView, DistributedView,
//! and Vector, plus the shape/contiguity errors each is expected to throw.

#include "UnitTest++.h"

#include <string>
#include <vector>

#include "AmanziComm.hh"
#include "AmanziTypes.hh"
#include "AmanziVector.hh"
#include "dbc.hh"
#include "errors.hh"
#include "FileHDF5.hh"

// FileHDF5's own error paths (Errors::Message, DBC::Assertion) are tested
// below by deliberately triggering HDF5-level failures (e.g. opening a
// dataset with a mismatched shape) -- HDF5's default error handler prints
// its internal error stack to stderr for every one of these, which is
// expected noise, not a sign anything is wrong.  CHECK_THROW_QUIET wraps
// CHECK_THROW in HDF5's own H5E_BEGIN_TRY/H5E_END_TRY (hdf5.h, included via
// FileHDF5.hh) to suppress just that expected printing for the wrapped
// expression; it saves and restores HDF5's current handler on every use, so
// an HDF5-level error anywhere else in the suite still prints normally.
#define CHECK_THROW_QUIET(expression, ExpectedExceptionType)                                       \
  do {                                                                                             \
    H5E_BEGIN_TRY                                                                                  \
    { CHECK_THROW(expression, ExpectedExceptionType); }                                            \
    H5E_END_TRY                                                                                    \
  } while (0)

SUITE(FILE_HDF5)
{
  using namespace Amanzi;

  //
  // Attribute -- scalar, at a group given by path (default root).
  //
  TEST(ATTR_WRITEREAD_SCALAR)
  {
    double sixone(6.1);
    int seven(7);
    bool truthy(true);
    long bignum(123456789L);
    std::string name("hello");

    {
      FileHDF5 out(getDefaultComm(), "test_hdf5_attr.h5", FILE_CREATE);
      out.writeAttribute("six_one", "/", sixone);
      out.writeAttribute("seven", "/", seven);
      out.writeAttribute("truthy", "/", truthy);
      out.writeAttribute("bignum", "/", bignum);
      out.writeAttribute("name", "/", name);
    }
    {
      FileHDF5 in(getDefaultComm(), "test_hdf5_attr.h5", FILE_READONLY);
      CHECK_CLOSE(sixone, in.readAttribute<double>("six_one"), 1.e-10);
      CHECK_EQUAL(seven, in.readAttribute<int>("seven"));
      CHECK_EQUAL(truthy, in.readAttribute<bool>("truthy"));
      CHECK_EQUAL(bignum, in.readAttribute<long>("bignum"));
      CHECK_EQUAL(name, in.readAttribute<std::string>("name"));
    }
  }

  // Attributes scoped to a non-root group, including a path encoded
  // directly in the attribute name (as OutputXDMF does for per-cycle
  // attributes, e.g. "42/time") -- both forms must create the group and
  // land the attribute in it.
  TEST(ATTR_WRITEREAD_GROUP_PATH)
  {
    auto comm = getDefaultComm();
    double time_42(4.2), time_43(4.3);

    {
      FileHDF5 out(comm, "test_hdf5_attr_group.h5", FILE_CREATE);
      out.writeAttribute("time", "42/", time_42);
      out.writeAttribute<double>("43/time", "/", time_43);
    }
    {
      FileHDF5 in(comm, "test_hdf5_attr_group.h5", FILE_READONLY);
      CHECK_CLOSE(time_42, in.readAttribute<double>("time", "42/"), 1.e-10);
      CHECK_CLOSE(time_43, in.readAttribute<double>("43/time"), 1.e-10);
    }
  }


  //
  // LocalView -- scalar, 1D, or 2D, identical on every rank.
  //
  TEST(LOCALVIEW_WRITEREAD_SCALAR)
  {
    auto comm = getDefaultComm();
    double dval(3.14);
    int ival(42);

    {
      FileHDF5 out(comm, "test_hdf5_localview_scalar.h5", FILE_CREATE);
      out.writeLocalView("dval", "/", dval);
      out.writeLocalView("ival", "/", ival);
    }
    {
      FileHDF5 in(comm, "test_hdf5_localview_scalar.h5", FILE_READONLY);
      CHECK_CLOSE(dval, in.readLocalView<double>("dval", "/"), 1.e-10);
      CHECK_EQUAL(ival, in.readLocalView<int>("ival", "/"));
    }
  }

  TEST(LOCALVIEW_WRITEREAD_1D)
  {
    auto comm = getDefaultComm();
    int n = 6;
    Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> vals("vals", n);
    for (int i = 0; i != n; ++i) vals(i) = 10.0 * i;

    {
      FileHDF5 out(comm, "test_hdf5_localview_1d.h5", FILE_CREATE);
      out.writeLocalView("myvar", "/", vals);
    }
    {
      FileHDF5 in(comm, "test_hdf5_localview_1d.h5", FILE_READONLY);
      Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> vals_in("vals_in", n);
      in.readLocalView("myvar", "/", vals_in);
      for (int i = 0; i != n; ++i) CHECK_CLOSE(vals(i), vals_in(i), 1.e-10);
    }
  }

  TEST(LOCALVIEW_WRITEREAD_2D)
  {
    auto comm = getDefaultComm();
    int n = 4, ndofs = 3;
    Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::HostSpace> vals("vals", n, ndofs);
    for (int i = 0; i != n; ++i)
      for (int j = 0; j != ndofs; ++j) vals(i, j) = 10 * i + j;

    {
      FileHDF5 out(comm, "test_hdf5_localview_2d.h5", FILE_CREATE);
      out.writeLocalView("myvar", "/", vals);
    }
    {
      FileHDF5 in(comm, "test_hdf5_localview_2d.h5", FILE_READONLY);
      Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::HostSpace> vals_in("vals_in", n, ndofs);
      in.readLocalView("myvar", "/", vals_in);
      for (int i = 0; i != n; ++i)
        for (int j = 0; j != ndofs; ++j) CHECK_CLOSE(vals(i, j), vals_in(i, j), 1.e-10);
    }
  }

  // A default-constructed (0-size) view passed to readLocalView must be
  // resized to the on-disk shape rather than throwing a shape-mismatch
  // error -- this is the mechanism callers use when they don't yet know
  // the shape of what they're about to read.
  TEST(LOCALVIEW_READ_1D_RESIZES_DEFAULT_CONSTRUCTED_VIEW)
  {
    auto comm = getDefaultComm();
    int n = 6;
    Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> vals("vals", n);
    for (int i = 0; i != n; ++i) vals(i) = 10.0 * i;

    {
      FileHDF5 out(comm, "test_hdf5_localview_1d_resize.h5", FILE_CREATE);
      out.writeLocalView("myvar", "/", vals);
    }
    {
      FileHDF5 in(comm, "test_hdf5_localview_1d_resize.h5", FILE_READONLY);
      Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> vals_in; // default-constructed
      CHECK_EQUAL(0, (int)vals_in.extent(0));
      in.readLocalView("myvar", "/", vals_in);
      CHECK_EQUAL(n, (int)vals_in.extent(0));
      for (int i = 0; i != n; ++i) CHECK_CLOSE(vals(i), vals_in(i), 1.e-10);
    }
  }

  TEST(LOCALVIEW_READ_2D_RESIZES_DEFAULT_CONSTRUCTED_VIEW)
  {
    auto comm = getDefaultComm();
    int n = 4, ndofs = 3;
    Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::HostSpace> vals("vals", n, ndofs);
    for (int i = 0; i != n; ++i)
      for (int j = 0; j != ndofs; ++j) vals(i, j) = 10 * i + j;

    {
      FileHDF5 out(comm, "test_hdf5_localview_2d_resize.h5", FILE_CREATE);
      out.writeLocalView("myvar", "/", vals);
    }
    {
      FileHDF5 in(comm, "test_hdf5_localview_2d_resize.h5", FILE_READONLY);
      Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::HostSpace> vals_in; // default-constructed
      CHECK_EQUAL(0, (int)vals_in.extent(0));
      CHECK_EQUAL(0, (int)vals_in.extent(1));
      in.readLocalView("myvar", "/", vals_in);
      CHECK_EQUAL(n, (int)vals_in.extent(0));
      CHECK_EQUAL(ndofs, (int)vals_in.extent(1));
      for (int i = 0; i != n; ++i)
        for (int j = 0; j != ndofs; ++j) CHECK_CLOSE(vals(i, j), vals_in(i, j), 1.e-10);
    }
  }

  // LocalView, paged by a single time_index -- on-disk array has one more
  // leading dimension than value.
  TEST(LOCALVIEW_WRITEREAD_TIME_INDEX_SCALAR)
  {
    auto comm = getDefaultComm();
    int ntimes = 4;

    {
      FileHDF5 out(comm, "test_hdf5_localview_time_index_scalar.h5", FILE_CREATE);
      for (std::size_t t = 0; t != (std::size_t)ntimes; ++t)
        out.writeLocalView("myvar", "/", 100.0 * t, t);
    }
    {
      FileHDF5 in(comm, "test_hdf5_localview_time_index_scalar.h5", FILE_READONLY);
      for (std::size_t t = 0; t != (std::size_t)ntimes; ++t)
        CHECK_CLOSE(100.0 * t, in.readLocalView<double>("myvar", "/", t), 1.e-10);
    }
  }

  TEST(LOCALVIEW_WRITEREAD_TIME_INDEX_1D)
  {
    auto comm = getDefaultComm();
    int n = 3, ntimes = 4;

    {
      FileHDF5 out(comm, "test_hdf5_localview_time_index_1d.h5", FILE_CREATE);
      for (int t = 0; t != ntimes; ++t) {
        Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> vals("vals", n);
        for (int i = 0; i != n; ++i) vals(i) = 1000 * t + i;
        out.writeLocalView("myvar", "/", vals, (std::size_t)t);
      }
    }
    {
      FileHDF5 in(comm, "test_hdf5_localview_time_index_1d.h5", FILE_READONLY);
      for (int t = 0; t != ntimes; ++t) {
        Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> vals_in("vals_in", n);
        in.readLocalView("myvar", "/", vals_in, (std::size_t)t);
        for (int i = 0; i != n; ++i) CHECK_CLOSE(1000 * t + i, vals_in(i), 1.e-10);
      }
    }
  }

  TEST(LOCALVIEW_READ_TIME_INDEX_1D_RESIZES_DEFAULT_CONSTRUCTED_VIEW)
  {
    auto comm = getDefaultComm();
    int n = 3, ntimes = 4;

    {
      FileHDF5 out(comm, "test_hdf5_localview_time_index_1d_resize.h5", FILE_CREATE);
      for (int t = 0; t != ntimes; ++t) {
        Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> vals("vals", n);
        for (int i = 0; i != n; ++i) vals(i) = 1000 * t + i;
        out.writeLocalView("myvar", "/", vals, (std::size_t)t);
      }
    }
    {
      FileHDF5 in(comm, "test_hdf5_localview_time_index_1d_resize.h5", FILE_READONLY);
      for (int t = 0; t != ntimes; ++t) {
        Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace>
          vals_in; // default-constructed
        CHECK_EQUAL(0, (int)vals_in.extent(0));
        in.readLocalView("myvar", "/", vals_in, (std::size_t)t);
        CHECK_EQUAL(n, (int)vals_in.extent(0));
        for (int i = 0; i != n; ++i) CHECK_CLOSE(1000 * t + i, vals_in(i), 1.e-10);
      }
    }
  }

  TEST(LOCALVIEW_WRITEREAD_TIME_INDEX_2D)
  {
    auto comm = getDefaultComm();
    int n = 3, ndofs = 2, ntimes = 3;

    {
      FileHDF5 out(comm, "test_hdf5_localview_time_index_2d.h5", FILE_CREATE);
      for (int t = 0; t != ntimes; ++t) {
        Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::HostSpace> vals("vals", n, ndofs);
        for (int i = 0; i != n; ++i)
          for (int j = 0; j != ndofs; ++j) vals(i, j) = 1000 * t + 10 * i + j;
        out.writeLocalView("myvar", "/", vals, (std::size_t)t);
      }
    }
    {
      FileHDF5 in(comm, "test_hdf5_localview_time_index_2d.h5", FILE_READONLY);
      for (int t = 0; t != ntimes; ++t) {
        Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::HostSpace> vals_in("vals_in", n, ndofs);
        in.readLocalView("myvar", "/", vals_in, (std::size_t)t);
        for (int i = 0; i != n; ++i)
          for (int j = 0; j != ndofs; ++j)
            CHECK_CLOSE(1000 * t + 10 * i + j, vals_in(i, j), 1.e-10);
      }
    }
  }

  // LocalView, paged by a time_chunk range -- on-disk array has the same
  // rank as value, value's 0th dimension being the page length.
  TEST(LOCALVIEW_WRITEREAD_TIME_CHUNK_1D)
  {
    auto comm = getDefaultComm();
    int n = 3, ntimes = 5;

    Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::HostSpace> buf("buf", ntimes, n);
    for (int t = 0; t != ntimes; ++t)
      for (int i = 0; i != n; ++i) buf(t, i) = 1000 * t + i;

    {
      FileHDF5 out(comm, "test_hdf5_localview_time_chunk_1d.h5", FILE_CREATE);
      out.writeLocalView("myvar", "/", buf, FileHDF5::TimeRange{ 0, (std::size_t)ntimes });
    }
    {
      FileHDF5 in(comm, "test_hdf5_localview_time_chunk_1d.h5", FILE_READONLY);
      Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::HostSpace> buf_in("buf_in", ntimes, n);
      in.readLocalView("myvar", "/", buf_in, FileHDF5::TimeRange{ 0, (std::size_t)ntimes });
      for (int t = 0; t != ntimes; ++t)
        for (int i = 0; i != n; ++i) CHECK_CLOSE(buf(t, i), buf_in(t, i), 1.e-10);
    }
  }

  TEST(LOCALVIEW_READ_TIME_CHUNK_1D_RESIZES_DEFAULT_CONSTRUCTED_VIEW)
  {
    auto comm = getDefaultComm();
    int n = 3, ntimes = 5;

    Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::HostSpace> buf("buf", ntimes, n);
    for (int t = 0; t != ntimes; ++t)
      for (int i = 0; i != n; ++i) buf(t, i) = 1000 * t + i;

    {
      FileHDF5 out(comm, "test_hdf5_localview_time_chunk_1d_resize.h5", FILE_CREATE);
      out.writeLocalView("myvar", "/", buf, FileHDF5::TimeRange{ 0, (std::size_t)ntimes });
    }
    {
      FileHDF5 in(comm, "test_hdf5_localview_time_chunk_1d_resize.h5", FILE_READONLY);
      Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::HostSpace> buf_in; // default-constructed
      CHECK_EQUAL(0, (int)buf_in.extent(0));
      CHECK_EQUAL(0, (int)buf_in.extent(1));
      in.readLocalView("myvar", "/", buf_in, FileHDF5::TimeRange{ 0, (std::size_t)ntimes });
      CHECK_EQUAL(ntimes, (int)buf_in.extent(0));
      CHECK_EQUAL(n, (int)buf_in.extent(1));
      for (int t = 0; t != ntimes; ++t)
        for (int i = 0; i != n; ++i) CHECK_CLOSE(buf(t, i), buf_in(t, i), 1.e-10);
    }
  }

  TEST(LOCALVIEW_WRITEREAD_TIME_CHUNK_2D)
  {
    auto comm = getDefaultComm();
    int n = 3, ndofs = 2, ntimes = 4;

    Kokkos::View<double***, Kokkos::LayoutRight, Kokkos::HostSpace> buf("buf", ntimes, n, ndofs);
    for (int t = 0; t != ntimes; ++t)
      for (int i = 0; i != n; ++i)
        for (int j = 0; j != ndofs; ++j) buf(t, i, j) = 1000 * t + 10 * i + j;

    {
      FileHDF5 out(comm, "test_hdf5_localview_time_chunk_2d.h5", FILE_CREATE);
      out.writeLocalView("myvar", "/", buf, FileHDF5::TimeRange{ 0, (std::size_t)ntimes });
    }
    {
      FileHDF5 in(comm, "test_hdf5_localview_time_chunk_2d.h5", FILE_READONLY);
      Kokkos::View<double***, Kokkos::LayoutRight, Kokkos::HostSpace> buf_in(
        "buf_in", ntimes, n, ndofs);
      in.readLocalView("myvar", "/", buf_in, FileHDF5::TimeRange{ 0, (std::size_t)ntimes });
      for (int t = 0; t != ntimes; ++t)
        for (int i = 0; i != n; ++i)
          for (int j = 0; j != ndofs; ++j) CHECK_CLOSE(buf(t, i, j), buf_in(t, i, j), 1.e-10);
    }
  }

  //
  // Open-ended appending: a dataset's time dimension must support being
  // extended indefinitely, in pieces, not just written once in full, and
  // the two granularities (time_index, time_chunk) must interoperate.
  //
  TEST(LOCALVIEW_MULTI_CHUNK_APPEND)
  {
    auto comm = getDefaultComm();
    int n = 3, ntimes = 7;

    Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::HostSpace> buf("buf", ntimes, n);
    for (int t = 0; t != ntimes; ++t)
      for (int i = 0; i != n; ++i) buf(t, i) = 1000 * t + i;

    {
      FileHDF5 out(comm, "test_hdf5_localview_multi_chunk_append.h5", FILE_CREATE);
      std::vector<FileHDF5::TimeRange> pieces{ { 0, 2 }, { 2, 5 }, { 5, 7 } };
      for (const auto& piece : pieces) {
        int count = piece.second - piece.first;
        Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::HostSpace> slice("slice", count, n);
        for (int t = 0; t != count; ++t)
          for (int i = 0; i != n; ++i) slice(t, i) = buf(piece.first + t, i);
        out.writeLocalView("myvar", "/", slice, piece);
      }
    }

    Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::HostSpace> buf_in("buf_in", ntimes, n);
    {
      FileHDF5 in(comm, "test_hdf5_localview_multi_chunk_append.h5", FILE_READONLY);
      in.readLocalView("myvar", "/", buf_in, FileHDF5::TimeRange{ 0, (std::size_t)ntimes });
    }
    for (int t = 0; t != ntimes; ++t)
      for (int i = 0; i != n; ++i) CHECK_CLOSE(buf(t, i), buf_in(t, i), 1.e-10);
  }

  TEST(LOCALVIEW_WRITE_CHUNK_READ_INDEX)
  {
    auto comm = getDefaultComm();
    int n = 3, ntimes = 4;

    Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::HostSpace> buf("buf", ntimes, n);
    for (int t = 0; t != ntimes; ++t)
      for (int i = 0; i != n; ++i) buf(t, i) = 1000 * t + i;

    {
      FileHDF5 out(comm, "test_hdf5_localview_write_chunk_read_index.h5", FILE_CREATE);
      out.writeLocalView("myvar", "/", buf, FileHDF5::TimeRange{ 0, (std::size_t)ntimes });
    }
    {
      FileHDF5 in(comm, "test_hdf5_localview_write_chunk_read_index.h5", FILE_READONLY);
      for (int t = 0; t != ntimes; ++t) {
        Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> s_in("s_in", n);
        in.readLocalView("myvar", "/", s_in, (std::size_t)t);
        for (int i = 0; i != n; ++i) CHECK_CLOSE(buf(t, i), s_in(i), 1.e-10);
      }
    }
  }

  TEST(LOCALVIEW_WRITE_INDEX_READ_CHUNK)
  {
    auto comm = getDefaultComm();
    int n = 3, ntimes = 4;

    Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::HostSpace> buf("buf", ntimes, n);
    for (int t = 0; t != ntimes; ++t)
      for (int i = 0; i != n; ++i) buf(t, i) = 1000 * t + i;

    {
      FileHDF5 out(comm, "test_hdf5_localview_write_index_read_chunk.h5", FILE_CREATE);
      for (int t = 0; t != ntimes; ++t) {
        Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> slice("slice", n);
        for (int i = 0; i != n; ++i) slice(i) = buf(t, i);
        out.writeLocalView("myvar", "/", slice, (std::size_t)t);
      }
    }

    Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::HostSpace> buf_in("buf_in", ntimes, n);
    {
      FileHDF5 in(comm, "test_hdf5_localview_write_index_read_chunk.h5", FILE_READONLY);
      in.readLocalView("myvar", "/", buf_in, FileHDF5::TimeRange{ 0, (std::size_t)ntimes });
    }
    for (int t = 0; t != ntimes; ++t)
      for (int i = 0; i != n; ++i) CHECK_CLOSE(buf(t, i), buf_in(t, i), 1.e-10);
  }

  TEST(LOCALVIEW_REOPEN_AND_APPEND)
  {
    // Simulates indefinite appending across separate runs: each "run" opens
    // a fresh FileHDF5 object in FILE_READWRITE mode (constructor opens,
    // destructor closes), writes one more time index past the previous
    // run's end, then closes.
    auto comm = getDefaultComm();
    int n = 3, ntimes = 5;

    Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::HostSpace> buf("buf", ntimes, n);
    for (int t = 0; t != ntimes; ++t)
      for (int i = 0; i != n; ++i) buf(t, i) = 1000 * t + i;

    {
      FileHDF5 create(comm, "test_hdf5_localview_reopen_append.h5", FILE_CREATE);
    }
    for (int t = 0; t != ntimes; ++t) {
      FileHDF5 out(comm, "test_hdf5_localview_reopen_append.h5", FILE_READWRITE);
      Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> slice("slice", n);
      for (int i = 0; i != n; ++i) slice(i) = buf(t, i);
      out.writeLocalView("myvar", "/", slice, (std::size_t)t);
    }

    {
      FileHDF5 in(comm, "test_hdf5_localview_reopen_append.h5", FILE_READONLY);
      for (int t = 0; t != ntimes; ++t) {
        Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> s_in("s_in", n);
        in.readLocalView("myvar", "/", s_in, (std::size_t)t);
        for (int i = 0; i != n; ++i) CHECK_CLOSE(buf(t, i), s_in(i), 1.e-10);
      }
    }
  }


  //
  // DistributedView -- 1D, decomposed across ranks via ASCEMIO.
  //
  TEST(DISTRIBUTEDVIEW_WRITEREAD_1D)
  {
    auto comm = getDefaultComm();
    int rank = comm->getRank();
    int nlocal = 5;

    Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> vals("vals", nlocal);
    for (int i = 0; i != nlocal; ++i) vals(i) = 10 * rank + i;

    {
      FileHDF5 out(comm, "test_hdf5_distview_1d.h5", FILE_CREATE);
      out.writeDistributedView("myvar", "/", vals);
    }
    {
      FileHDF5 in(comm, "test_hdf5_distview_1d.h5", FILE_READONLY);
      Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> vals_in("vals_in", nlocal);
      in.readDistributedView("myvar", "/", vals_in);
      for (int i = 0; i != nlocal; ++i) CHECK_CLOSE(vals(i), vals_in(i), 1.e-10);
    }
  }

  TEST(DISTRIBUTEDVIEW_WRITEREAD_1D_INT)
  {
    auto comm = getDefaultComm();
    int rank = comm->getRank();
    int nlocal = 4;

    Kokkos::View<int*, Kokkos::LayoutRight, Kokkos::HostSpace> vals("vals", nlocal);
    for (int i = 0; i != nlocal; ++i) vals(i) = 100 * rank + i;

    {
      FileHDF5 out(comm, "test_hdf5_distview_1d_int.h5", FILE_CREATE);
      out.writeDistributedView("myvar", "/", vals);
    }
    {
      FileHDF5 in(comm, "test_hdf5_distview_1d_int.h5", FILE_READONLY);
      Kokkos::View<int*, Kokkos::LayoutRight, Kokkos::HostSpace> vals_in("vals_in", nlocal);
      in.readDistributedView("myvar", "/", vals_in);
      for (int i = 0; i != nlocal; ++i) CHECK_EQUAL(vals(i), vals_in(i));
    }
  }

  // DistributedView -- 2D [N, NDOFS], decomposed along dimension 0.
  TEST(DISTRIBUTEDVIEW_WRITEREAD_2D)
  {
    auto comm = getDefaultComm();
    int rank = comm->getRank();
    int nlocal = 4, ndofs = 3;

    Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::HostSpace> vals("vals", nlocal, ndofs);
    for (int i = 0; i != nlocal; ++i)
      for (int j = 0; j != ndofs; ++j) vals(i, j) = 100 * rank + 10 * i + j;

    {
      FileHDF5 out(comm, "test_hdf5_distview_2d.h5", FILE_CREATE);
      out.writeDistributedView("myvar", "/", vals);
    }
    {
      FileHDF5 in(comm, "test_hdf5_distview_2d.h5", FILE_READONLY);
      Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::HostSpace> vals_in(
        "vals_in", nlocal, ndofs);
      in.readDistributedView("myvar", "/", vals_in);
      for (int i = 0; i != nlocal; ++i)
        for (int j = 0; j != ndofs; ++j) CHECK_CLOSE(vals(i, j), vals_in(i, j), 1.e-10);
    }
  }


  //
  // Vector -- delegates to the 1D DistributedView.
  //
  TEST(VECTOR_WRITEREAD)
  {
    auto comm = getDefaultComm();
    int size = comm->getSize();
    int rank = comm->getRank();
    Map_ptr_type map = Teuchos::rcp(new Map_type(size * 3, 3, 0, comm));

    Vector_type vec(map);
    {
      auto vecv = vec.getLocalViewHost(Tpetra::Access::ReadWrite);
      for (int i = 0; i != (int)vecv.extent(0); ++i) vecv(i, 0) = i + rank * 3;
    }

    {
      FileHDF5 out(comm, "test_hdf5_vector.h5", FILE_CREATE);
      out.writeVector("myvec", "/", vec);
    }

    Vector_type vec_in(map);
    {
      FileHDF5 in(comm, "test_hdf5_vector.h5", FILE_READONLY);
      in.readVector("myvec", "/", vec_in);
    }

    vec_in.update(-1, vec, 1);
    CHECK_CLOSE(0.0, vec_in.norm2(), 1.e-10);
  }

  TEST(INTVECTOR_WRITEREAD)
  {
    auto comm = getDefaultComm();
    int size = comm->getSize();
    int rank = comm->getRank();
    Map_ptr_type map = Teuchos::rcp(new Map_type(size * 4, 4, 0, comm));

    IntVector_type vec(map);
    {
      auto vecv = vec.getLocalViewHost(Tpetra::Access::ReadWrite);
      for (int i = 0; i != (int)vecv.extent(0); ++i) vecv(i, 0) = i + rank * 4;
    }

    {
      FileHDF5 out(comm, "test_hdf5_intvector.h5", FILE_CREATE);
      out.writeVector("myvec", "/", vec);
    }

    IntVector_type vec_in(map);
    {
      FileHDF5 in(comm, "test_hdf5_intvector.h5", FILE_READONLY);
      in.readVector("myvec", "/", vec_in);
    }

    auto vecv = vec.getLocalViewHost(Tpetra::Access::ReadOnly);
    auto vecv_in = vec_in.getLocalViewHost(Tpetra::Access::ReadOnly);
    for (int i = 0; i != (int)vecv.extent(0); ++i) CHECK_EQUAL(vecv(i, 0), vecv_in(i, 0));
  }

  // readVector only needs a Vector_type-shaped argument -- one column of a
  // MultiVector (via getVectorNonConst) is exactly that, so a vector written
  // standalone must be readable directly into a single dof of a multi-dof
  // MultiVector, without disturbing that MultiVector's other columns.
  TEST(VECTOR_READ_INTO_MULTIVECTOR_COLUMN)
  {
    auto comm = getDefaultComm();
    int size = comm->getSize();
    int rank = comm->getRank();
    Map_ptr_type map = Teuchos::rcp(new Map_type(size * 3, 3, 0, comm));

    Vector_type vec(map);
    {
      auto vecv = vec.getLocalViewHost(Tpetra::Access::ReadWrite);
      for (int i = 0; i != (int)vecv.extent(0); ++i) vecv(i, 0) = i + rank * 3;
    }

    {
      FileHDF5 out(comm, "test_hdf5_vector_for_multivector.h5", FILE_CREATE);
      out.writeVector("myvec", "/", vec);
    }

    int ndofs = 4;
    int target_dof = 2;
    MultiVector_type mvec(map, ndofs);
    {
      // fill every column with an easily-recognized, distinct pattern, so
      // reading into target_dof must not disturb the others.
      auto mvecv = mvec.getLocalViewHost(Tpetra::Access::ReadWrite);
      for (int i = 0; i != (int)mvecv.extent(0); ++i)
        for (int j = 0; j != ndofs; ++j) mvecv(i, j) = -1000 * (j + 1) - i - rank * 3;
    }

    {
      FileHDF5 in(comm, "test_hdf5_vector_for_multivector.h5", FILE_READONLY);
      Teuchos::RCP<Vector_type> col = mvec.getVectorNonConst(target_dof);
      in.readVector("myvec", "/", *col);
    }

    auto vecv = vec.getLocalViewHost(Tpetra::Access::ReadOnly);
    auto mvecv = mvec.getLocalViewHost(Tpetra::Access::ReadOnly);
    for (int i = 0; i != (int)vecv.extent(0); ++i) {
      // the targeted column now matches the standalone vector exactly.
      CHECK_CLOSE(vecv(i, 0), mvecv(i, target_dof), 1.e-10);
      // the untouched columns still hold their original fill pattern.
      for (int j = 0; j != ndofs; ++j) {
        if (j == target_dof) continue;
        CHECK_CLOSE(-1000 * (j + 1) - i - rank * 3, mvecv(i, j), 1.e-10);
      }
    }
  }


  //
  // Error paths -- incorrectly shaped on-disk data must throw
  // Errors::Message when reading.
  //
  TEST(READ_LOCALVIEW_1D_WRONG_LENGTH_THROWS)
  {
    auto comm = getDefaultComm();
    Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> vals("vals", 5);
    for (int i = 0; i != 5; ++i) vals(i) = i;

    {
      FileHDF5 out(comm, "test_hdf5_localview_1d_wrong_length.h5", FILE_CREATE);
      out.writeLocalView("myvar", "/", vals);
    }
    {
      FileHDF5 in(comm, "test_hdf5_localview_1d_wrong_length.h5", FILE_READONLY);
      Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> vals_in(
        "vals_in", 4); // on-disk is length 5
      CHECK_THROW_QUIET(in.readLocalView("myvar", "/", vals_in), Errors::Message);
    }
  }

  TEST(READ_LOCALVIEW_RANK_MISMATCH_THROWS)
  {
    auto comm = getDefaultComm();
    Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> vals("vals", 6);
    for (int i = 0; i != 6; ++i) vals(i) = i;

    {
      FileHDF5 out(comm, "test_hdf5_localview_rank_mismatch.h5", FILE_CREATE);
      out.writeLocalView("myvar", "/", vals); // on-disk rank 1
    }
    {
      FileHDF5 in(comm, "test_hdf5_localview_rank_mismatch.h5", FILE_READONLY);
      Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::HostSpace> vals_in(
        "vals_in", 2, 3); // same span, rank 2
      CHECK_THROW_QUIET(in.readLocalView("myvar", "/", vals_in), Errors::Message);
    }
  }

  TEST(READ_LOCALVIEW_2D_WRONG_SHAPE_THROWS)
  {
    auto comm = getDefaultComm();
    Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::HostSpace> vals("vals", 4, 3);
    for (int i = 0; i != 4; ++i)
      for (int j = 0; j != 3; ++j) vals(i, j) = 10 * i + j;

    {
      FileHDF5 out(comm, "test_hdf5_localview_2d_wrong_shape.h5", FILE_CREATE);
      out.writeLocalView("myvar", "/", vals); // on-disk [4, 3]
    }
    {
      FileHDF5 in(comm, "test_hdf5_localview_2d_wrong_shape.h5", FILE_READONLY);
      Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::HostSpace> vals_in(
        "vals_in", 4, 2); // dim 1 mismatch
      CHECK_THROW_QUIET(in.readLocalView("myvar", "/", vals_in), Errors::Message);
    }
  }

  TEST(READ_LOCALVIEW_TIME_INDEX_WRONG_SHAPE_THROWS)
  {
    auto comm = getDefaultComm();
    Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> vals("vals", 3);
    for (int i = 0; i != 3; ++i) vals(i) = i;

    {
      FileHDF5 out(comm, "test_hdf5_localview_time_index_wrong_shape.h5", FILE_CREATE);
      out.writeLocalView("myvar", "/", vals, (std::size_t)0); // on-disk [1, 3]
    }
    {
      FileHDF5 in(comm, "test_hdf5_localview_time_index_wrong_shape.h5", FILE_READONLY);
      Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> vals_in("vals_in",
                                                                            4); // dim 1 mismatch
      CHECK_THROW_QUIET(in.readLocalView("myvar", "/", vals_in, (std::size_t)0), Errors::Message);
    }
  }

  TEST(READ_LOCALVIEW_TIME_CHUNK_WRONG_SHAPE_THROWS)
  {
    auto comm = getDefaultComm();
    Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::HostSpace> buf("buf", 4, 3);
    for (int t = 0; t != 4; ++t)
      for (int i = 0; i != 3; ++i) buf(t, i) = 10 * t + i;

    {
      FileHDF5 out(comm, "test_hdf5_localview_time_chunk_wrong_shape.h5", FILE_CREATE);
      out.writeLocalView("myvar", "/", buf, FileHDF5::TimeRange{ 0, 4 }); // on-disk [4, 3]
    }
    {
      FileHDF5 in(comm, "test_hdf5_localview_time_chunk_wrong_shape.h5", FILE_READONLY);
      Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::HostSpace> buf_in(
        "buf_in", 4, 2); // dim 1 mismatch
      CHECK_THROW_QUIET(in.readLocalView("myvar", "/", buf_in, (FileHDF5::TimeRange{ 0, 4 })),
                        Errors::Message);
    }
  }

  TEST(READ_DISTRIBUTEDVIEW_1D_WRONG_GLOBAL_LENGTH_THROWS)
  {
    auto comm = getDefaultComm();
    int rank = comm->getRank();
    int nlocal = 5;
    Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> vals("vals", nlocal);
    for (int i = 0; i != nlocal; ++i) vals(i) = 10 * rank + i;

    {
      FileHDF5 out(comm, "test_hdf5_distview_1d_wrong_length.h5", FILE_CREATE);
      out.writeDistributedView("myvar", "/", vals);
    }
    {
      FileHDF5 in(comm, "test_hdf5_distview_1d_wrong_length.h5", FILE_READONLY);
      Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> vals_in(
        "vals_in", nlocal + 1); // wrong global length
      CHECK_THROW_QUIET(in.readDistributedView("myvar", "/", vals_in), Errors::Message);
    }
  }

  TEST(READ_DISTRIBUTEDVIEW_2D_WRONG_SHAPE_THROWS)
  {
    auto comm = getDefaultComm();
    int rank = comm->getRank();
    int nlocal = 4, ndofs = 3;
    Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::HostSpace> vals("vals", nlocal, ndofs);
    for (int i = 0; i != nlocal; ++i)
      for (int j = 0; j != ndofs; ++j) vals(i, j) = 100 * rank + 10 * i + j;

    {
      FileHDF5 out(comm, "test_hdf5_distview_2d_wrong_shape.h5", FILE_CREATE);
      out.writeDistributedView("myvar", "/", vals); // on-disk [*, 3]
    }
    {
      FileHDF5 in(comm, "test_hdf5_distview_2d_wrong_shape.h5", FILE_READONLY);
      Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::HostSpace> vals_in(
        "vals_in", nlocal, ndofs + 1); // dim 1 mismatch
      CHECK_THROW_QUIET(in.readDistributedView("myvar", "/", vals_in), Errors::Message);
    }
  }

  // A 1D DistributedView writes a rank-2 on-disk dataset ([global_length, 1]
  // -- ASCEMIO's convention), so reading it back with a 2D view requires
  // the matching trailing dimension; a mismatched NDOFS must still throw.
  TEST(READ_DISTRIBUTEDVIEW_2D_FROM_1D_WRONG_NDOFS_THROWS)
  {
    auto comm = getDefaultComm();
    int rank = comm->getRank();
    int nlocal = 5;
    Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> vals("vals", nlocal);
    for (int i = 0; i != nlocal; ++i) vals(i) = 10 * rank + i;

    {
      FileHDF5 out(comm, "test_hdf5_distview_1d_as_2d.h5", FILE_CREATE);
      out.writeDistributedView("myvar", "/", vals); // on-disk [*, 1]
    }
    {
      FileHDF5 in(comm, "test_hdf5_distview_1d_as_2d.h5", FILE_READONLY);
      Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::HostSpace> vals_in(
        "vals_in", nlocal, 2); // NDOFS mismatch (on-disk NDOFS is 1)
      CHECK_THROW_QUIET(in.readDistributedView("myvar", "/", vals_in), Errors::Message);
    }
  }

  // Transposed shape: on-disk is [N, NDOFS] with N != NDOFS; reading with
  // the dimensions swapped has the same total element count but must still
  // throw, since the two dimensions are checked independently.
  TEST(READ_LOCALVIEW_2D_TRANSPOSED_SHAPE_THROWS)
  {
    auto comm = getDefaultComm();
    int n = 5, ndofs = 2;
    Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::HostSpace> vals("vals", n, ndofs);
    for (int i = 0; i != n; ++i)
      for (int j = 0; j != ndofs; ++j) vals(i, j) = 10 * i + j;

    {
      FileHDF5 out(comm, "test_hdf5_localview_2d_transposed.h5", FILE_CREATE);
      out.writeLocalView("myvar", "/", vals); // on-disk [5, 2]
    }
    {
      FileHDF5 in(comm, "test_hdf5_localview_2d_transposed.h5", FILE_READONLY);
      Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::HostSpace> vals_in(
        "vals_in", ndofs, n); // [2, 5] -- transposed, same span
      CHECK_THROW_QUIET(in.readLocalView("myvar", "/", vals_in), Errors::Message);
    }
  }

  TEST(READ_DISTRIBUTEDVIEW_2D_TRANSPOSED_SHAPE_THROWS)
  {
    auto comm = getDefaultComm();
    int rank = comm->getRank();
    int nlocal = 5, ndofs = 2;
    Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::HostSpace> vals("vals", nlocal, ndofs);
    for (int i = 0; i != nlocal; ++i)
      for (int j = 0; j != ndofs; ++j) vals(i, j) = 100 * rank + 10 * i + j;

    {
      FileHDF5 out(comm, "test_hdf5_distview_2d_transposed.h5", FILE_CREATE);
      out.writeDistributedView("myvar", "/", vals); // on-disk [*, 2], * == nlocal * nranks
    }
    {
      FileHDF5 in(comm, "test_hdf5_distview_2d_transposed.h5", FILE_READONLY);
      // same local element count (nlocal * ndofs), dimensions swapped.
      Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::HostSpace> vals_in(
        "vals_in", ndofs, nlocal);
      CHECK_THROW_QUIET(in.readDistributedView("myvar", "/", vals_in), Errors::Message);
    }
  }

  // Too few / too many entries: reading a shorter or longer 1D buffer than
  // what was written must throw in both directions, not just one.
  TEST(READ_LOCALVIEW_1D_TOO_FEW_ENTRIES_THROWS)
  {
    auto comm = getDefaultComm();
    Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> vals("vals", 8);
    for (int i = 0; i != 8; ++i) vals(i) = i;

    {
      FileHDF5 out(comm, "test_hdf5_localview_1d_too_few.h5", FILE_CREATE);
      out.writeLocalView("myvar", "/", vals); // on-disk length 8
    }
    {
      FileHDF5 in(comm, "test_hdf5_localview_1d_too_few.h5", FILE_READONLY);
      Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> vals_in("vals_in", 3);
      CHECK_THROW_QUIET(in.readLocalView("myvar", "/", vals_in), Errors::Message);
    }
  }

  TEST(READ_LOCALVIEW_1D_TOO_MANY_ENTRIES_THROWS)
  {
    auto comm = getDefaultComm();
    Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> vals("vals", 8);
    for (int i = 0; i != 8; ++i) vals(i) = i;

    {
      FileHDF5 out(comm, "test_hdf5_localview_1d_too_many.h5", FILE_CREATE);
      out.writeLocalView("myvar", "/", vals); // on-disk length 8
    }
    {
      FileHDF5 in(comm, "test_hdf5_localview_1d_too_many.h5", FILE_READONLY);
      Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> vals_in("vals_in", 20);
      CHECK_THROW_QUIET(in.readLocalView("myvar", "/", vals_in), Errors::Message);
    }
  }

  TEST(READ_DISTRIBUTEDVIEW_1D_TOO_FEW_ENTRIES_THROWS)
  {
    auto comm = getDefaultComm();
    int rank = comm->getRank();
    int nlocal = 6;
    Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> vals("vals", nlocal);
    for (int i = 0; i != nlocal; ++i) vals(i) = 10 * rank + i;

    {
      FileHDF5 out(comm, "test_hdf5_distview_1d_too_few.h5", FILE_CREATE);
      out.writeDistributedView("myvar", "/", vals);
    }
    {
      FileHDF5 in(comm, "test_hdf5_distview_1d_too_few.h5", FILE_READONLY);
      Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> vals_in(
        "vals_in", nlocal - 2); // reading fewer entries than were written
      CHECK_THROW_QUIET(in.readDistributedView("myvar", "/", vals_in), Errors::Message);
    }
  }

  TEST(READ_DISTRIBUTEDVIEW_1D_TOO_MANY_ENTRIES_THROWS)
  {
    auto comm = getDefaultComm();
    int rank = comm->getRank();
    int nlocal = 6;
    Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> vals("vals", nlocal);
    for (int i = 0; i != nlocal; ++i) vals(i) = 10 * rank + i;

    {
      FileHDF5 out(comm, "test_hdf5_distview_1d_too_many.h5", FILE_CREATE);
      out.writeDistributedView("myvar", "/", vals);
    }
    {
      FileHDF5 in(comm, "test_hdf5_distview_1d_too_many.h5", FILE_READONLY);
      Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> vals_in(
        "vals_in", nlocal + 3); // reading more entries than were written
      CHECK_THROW_QUIET(in.readDistributedView("myvar", "/", vals_in), Errors::Message);
    }
  }


  //
  // Error paths -- non-contiguous views must throw (AMANZI_ASSERT, i.e.
  // DBC::Assertion) rather than silently reading/writing a wrong buffer.
  //
  // Kokkos::subview() of a LayoutRight view "demotes" the result to
  // LayoutStride whenever it isn't trivially contiguous, so it cannot be
  // used to build a non-contiguous *LayoutRight* view -- it would simply
  // fail to compile against writeLocalView/writeDistributedView<2D>'s
  // LayoutRight constraint.  Instead, construct the LayoutRight view
  // directly with a padded (non-default) stride: Kokkos::LayoutRight's
  // stride pads the fastest-varying (last) dimension, so dimension 0 is
  // separated by more than extent(1) elements, making span_is_contiguous()
  // false while array_layout stays Kokkos::LayoutRight.
  TEST(WRITE_LOCALVIEW_2D_NONCONTIGUOUS_THROWS)
  {
    auto comm = getDefaultComm();
    int n = 3, ndofs = 2;
    Kokkos::LayoutRight padded_layout(n, ndofs);
    padded_layout.stride = ndofs + 2; // pad past the true row length
    Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::HostSpace> vals("vals", padded_layout);
    CHECK(!vals.span_is_contiguous());

    FileHDF5 out(comm, "test_hdf5_localview_2d_noncontiguous.h5", FILE_CREATE);
    CHECK_THROW_QUIET(out.writeLocalView("myvar", "/", vals), DBC::Assertion);
  }

  TEST(READ_LOCALVIEW_2D_NONCONTIGUOUS_THROWS)
  {
    auto comm = getDefaultComm();
    int n = 3, ndofs = 2;
    Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::HostSpace> vals("vals", n, ndofs);
    for (int i = 0; i != n; ++i)
      for (int j = 0; j != ndofs; ++j) vals(i, j) = 10 * i + j;

    {
      FileHDF5 out(comm, "test_hdf5_localview_2d_noncontiguous_read.h5", FILE_CREATE);
      out.writeLocalView("myvar", "/", vals);
    }

    Kokkos::LayoutRight padded_layout(n, ndofs);
    padded_layout.stride = ndofs + 2;
    Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::HostSpace> vals_in("vals_in",
                                                                           padded_layout);
    CHECK(!vals_in.span_is_contiguous());

    FileHDF5 in(comm, "test_hdf5_localview_2d_noncontiguous_read.h5", FILE_READONLY);
    CHECK_THROW_QUIET(in.readLocalView("myvar", "/", vals_in), DBC::Assertion);
  }

  TEST(WRITE_DISTRIBUTEDVIEW_2D_NONCONTIGUOUS_THROWS)
  {
    auto comm = getDefaultComm();
    int nlocal = 4, ndofs = 3;
    Kokkos::LayoutRight padded_layout(nlocal, ndofs);
    padded_layout.stride = ndofs + 1;
    Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::HostSpace> vals("vals", padded_layout);
    CHECK(!vals.span_is_contiguous());

    FileHDF5 out(comm, "test_hdf5_distview_2d_noncontiguous.h5", FILE_CREATE);
    CHECK_THROW_QUIET(out.writeDistributedView("myvar", "/", vals), DBC::Assertion);
  }

  TEST(WRITE_LOCALVIEW_TIME_CHUNK_NONCONTIGUOUS_THROWS)
  {
    auto comm = getDefaultComm();
    int n = 3, ndofs = 2, ntimes = 2;
    Kokkos::LayoutRight padded_layout(ntimes, n, ndofs);
    padded_layout.stride = n * ndofs + 2;
    Kokkos::View<double***, Kokkos::LayoutRight, Kokkos::HostSpace> vals("vals", padded_layout);
    CHECK(!vals.span_is_contiguous());

    FileHDF5 out(comm, "test_hdf5_localview_time_chunk_noncontiguous.h5", FILE_CREATE);
    CHECK_THROW_QUIET(
      out.writeLocalView("myvar", "/", vals, FileHDF5::TimeRange{ 0, (std::size_t)ntimes }),
      DBC::Assertion);
  }

  // DistributedView's rank-1 overload has no layout constraint (any layout
  // is as contiguous as any other for rank 1 in general), so a
  // non-contiguous 1D view needs an explicit LayoutStride, not LayoutRight
  // padding (LayoutRight is always contiguous at rank 1).
  TEST(WRITE_DISTRIBUTEDVIEW_1D_NONCONTIGUOUS_THROWS)
  {
    auto comm = getDefaultComm();
    int nlocal = 5;
    Kokkos::View<double*, Kokkos::LayoutStride, Kokkos::HostSpace> vals(
      "vals", Kokkos::LayoutStride(nlocal, 2));
    CHECK(!vals.span_is_contiguous());

    FileHDF5 out(comm, "test_hdf5_distview_1d_noncontiguous.h5", FILE_CREATE);
    CHECK_THROW_QUIET(out.writeDistributedView("myvar", "/", vals), DBC::Assertion);
  }

  TEST(READ_DISTRIBUTEDVIEW_1D_NONCONTIGUOUS_THROWS)
  {
    auto comm = getDefaultComm();
    int rank = comm->getRank();
    int nlocal = 5;
    Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> vals("vals", nlocal);
    for (int i = 0; i != nlocal; ++i) vals(i) = 10 * rank + i;

    {
      FileHDF5 out(comm, "test_hdf5_distview_1d_noncontiguous_read.h5", FILE_CREATE);
      out.writeDistributedView("myvar", "/", vals);
    }

    Kokkos::View<double*, Kokkos::LayoutStride, Kokkos::HostSpace> vals_in(
      "vals_in", Kokkos::LayoutStride(nlocal, 2));
    CHECK(!vals_in.span_is_contiguous());

    FileHDF5 in(comm, "test_hdf5_distview_1d_noncontiguous_read.h5", FILE_READONLY);
    CHECK_THROW_QUIET(in.readDistributedView("myvar", "/", vals_in), DBC::Assertion);
  }
}
