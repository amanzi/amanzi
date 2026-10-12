/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

//! Tests FileNetCDF's read/write API (attributes and Views).
/*
  Covers both of FileNetCDF.hh's shape families -- "packed" (decomposed_axis
  is the view's last axis) and "columnar" (a real NCELLS_PER_COLUMN axis
  kept separate from NCOLUMNS, decomposed_axis names NCOLUMNS's position,
  never the last axis).  There is no stored decomposition state on the
  object: each call's decomposed_axis, together with every rank's local
  extent at that axis, is all FileNetCDF needs (global length is a sum
  across ranks, offset an exclusive prefix sum).
*/

#include "UnitTest++.h"

#include <string>
#include <vector>

#include "AmanziComm.hh"
#include "AmanziTypes.hh"
#include "AmanziVector.hh"
#include "dbc.hh"
#include "errors.hh"
#include "FileNetCDF.hh"

extern "C"
{
#include "netcdf.h"
}

SUITE(FILE_NETCDF)
{
  using namespace Amanzi;

  TEST(ATTR_WRITEREAD_SCALAR)
  {
    double sixone(6.1);
    int seven(7);
    bool truthy(true);
    std::string name("hello");

    {
      FileNetCDF out(getDefaultComm(), "test_netcdf_attr.nc", FILE_CREATE);
      out.writeAttribute("six_one", sixone);
      out.writeAttribute("seven", seven);
      out.writeAttribute("truthy", truthy);
      out.writeAttribute("name", name);
    }
    {
      FileNetCDF in(getDefaultComm(), "test_netcdf_attr.nc", FILE_READONLY);
      CHECK_CLOSE(sixone, in.readAttribute<double>("six_one"), 1.e-10);
      CHECK_EQUAL(seven, in.readAttribute<int>("seven"));
      CHECK_EQUAL(truthy, in.readAttribute<bool>("truthy"));
      CHECK_EQUAL(name, in.readAttribute<std::string>("name"));
    }
  }

  TEST(ATTR_WRITEREAD_TIME_CHUNK)
  {
    // time-dependent attribute (e.g. "time"), written/read as a chunk.
    int ntimes = 6;
    Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> times("times", ntimes);
    for (int t = 0; t != ntimes; ++t) times(t) = 0.5 * t;

    {
      FileNetCDF out(getDefaultComm(), "test_netcdf_attr_time.nc", FILE_CREATE);
      out.writeAttribute("time", times, FileNetCDF::TimeRange{ 0, (std::size_t)ntimes });
    }
    {
      FileNetCDF in(getDefaultComm(), "test_netcdf_attr_time.nc", FILE_READONLY);
      Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> times_in("times_in", ntimes);
      in.readAttribute("time", times_in, FileNetCDF::TimeRange{ 0, (std::size_t)ntimes });
      for (int t = 0; t != ntimes; ++t) CHECK_CLOSE(times(t), times_in(t), 1.e-10);
    }
  }

  TEST(ATTR_WRITEREAD_TIME_INDEX)
  {
    // time-dependent attribute, written/read one index at a time.
    int ntimes = 4;

    {
      FileNetCDF out(getDefaultComm(), "test_netcdf_attr_time_index.nc", FILE_CREATE);
      for (int t = 0; t != ntimes; ++t) out.writeAttribute("time", 0.25 * t, t);
    }
    {
      FileNetCDF in(getDefaultComm(), "test_netcdf_attr_time_index.nc", FILE_READONLY);
      for (int t = 0; t != ntimes; ++t) {
        double val = in.readAttribute<double>("time", t);
        CHECK_CLOSE(0.25 * t, val, 1.e-10);
      }
    }
  }


  // Confirms a time-dependent attribute other than "time" itself shares
  // the "time" dimension on disk (not a dimension named after itself) --
  // every ATTR_WRITEREAD_TIME_* test above only ever writes under the name
  // "time", so this is the one case that actually exercises the dimension
  // name writeAttribute's TimeRange overload uses when defining a new
  // variable.
  TEST(ATTR_TIME_INDEX_SHARES_TIME_DIMENSION)
  {
    int ntimes = 3;
    std::string filename = "test_netcdf_attr_time_shared_dim.nc";

    {
      FileNetCDF out(getDefaultComm(), filename, FILE_CREATE);
      for (int t = 0; t != ntimes; ++t) {
        out.writeAttribute("time", 0.5 * t, t);
        out.writeAttribute("atmospheric_pressure", 101325.0 + t, t);
      }
    }

    int ncid;
    CHECK_EQUAL(NC_NOERR, nc_open(filename.c_str(), NC_NOWRITE, &ncid));

    // the "time" dimension exists (confirms it's literally named "time",
    // not just that the two variables below happen to share some dimension)
    int named_time_dimid;
    CHECK_EQUAL(NC_NOERR, nc_inq_dimid(ncid, "time", &named_time_dimid));

    int time_varid, atm_varid;
    CHECK_EQUAL(NC_NOERR, nc_inq_varid(ncid, "time", &time_varid));
    CHECK_EQUAL(NC_NOERR, nc_inq_varid(ncid, "atmospheric_pressure", &atm_varid));

    int time_dimid, atm_dimid;
    CHECK_EQUAL(NC_NOERR, nc_inq_vardimid(ncid, time_varid, &time_dimid));
    CHECK_EQUAL(NC_NOERR, nc_inq_vardimid(ncid, atm_varid, &atm_dimid));
    CHECK_EQUAL(time_dimid, atm_dimid);
    CHECK_EQUAL(named_time_dimid, atm_dimid);

    // and no separate "atmospheric_pressure" dimension was created
    int stray_dimid;
    CHECK_EQUAL(NC_EBADDIM, nc_inq_dimid(ncid, "atmospheric_pressure", &stray_dimid));

    nc_close(ncid);
  }

  //
  // Packed shapes: decomposed_axis is always the view's last axis.
  //
  TEST(VIEW_WRITEREAD_1D)
  {
    // shape [N]: Kokkos::View<T*>, decomposed_axis=0.
    auto comm = getDefaultComm();
    int rank = comm->getRank();
    int nranks = comm->getSize();
    int nlocal = 5;
    int offset = rank * nlocal;

    Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> vals("vals", nlocal);
    for (int i = 0; i != nlocal; ++i) vals(i) = 10 * (offset + i);

    {
      FileNetCDF out(comm, "test_netcdf_view_1d.nc", FILE_CREATE);
      out.writeDistributedView("myvar", vals, 0, { "n" });
    }
    {
      FileNetCDF in(comm, "test_netcdf_view_1d.nc", FILE_READONLY);
      Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> vals_in("vals_in", nlocal);
      in.readDistributedView("myvar", vals_in, 0, { "n" });
      for (int i = 0; i != nlocal; ++i) CHECK_CLOSE(vals(i), vals_in(i), 1.e-10);
    }
  }

  TEST(VIEW_WRITEREAD_PACKED_MULTIDOF)
  {
    // shape [NDOFS, N]: Kokkos::View<T**>, decomposed_axis=1.
    auto comm = getDefaultComm();
    int rank = comm->getRank();
    int nlocal = 4;
    int ndofs = 3;
    int offset = rank * nlocal;

    Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::HostSpace> vals("vals", ndofs, nlocal);
    for (int j = 0; j != ndofs; ++j)
      for (int i = 0; i != nlocal; ++i) vals(j, i) = 100 * (offset + i) + j;

    {
      FileNetCDF out(comm, "test_netcdf_view_packed_multidof.nc", FILE_CREATE);
      out.writeDistributedView("myvar", vals, 1, { "ndofs", "n" });
    }
    {
      FileNetCDF in(comm, "test_netcdf_view_packed_multidof.nc", FILE_READONLY);
      Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::HostSpace> vals_in(
        "vals_in", ndofs, nlocal);
      in.readDistributedView("myvar", vals_in, 1, { "ndofs", "n" });
      for (int j = 0; j != ndofs; ++j)
        for (int i = 0; i != nlocal; ++i) CHECK_CLOSE(vals(j, i), vals_in(j, i), 1.e-10);
    }
  }

  TEST(VIEW_WRITEREAD_PACKED_TIME_CHUNK_1DOF)
  {
    // shape [NTIMES, N]: Kokkos::View<T**>, decomposed_axis=1, with time_chunk.
    auto comm = getDefaultComm();
    int rank = comm->getRank();
    int nlocal = 3;
    int ntimes = 5;
    int offset = rank * nlocal;

    Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::HostSpace> vals("vals", ntimes, nlocal);
    for (int t = 0; t != ntimes; ++t)
      for (int i = 0; i != nlocal; ++i) vals(t, i) = 1000 * t + 10 * (offset + i);

    {
      FileNetCDF out(comm, "test_netcdf_view_packed_time_chunk.nc", FILE_CREATE);
      out.writeDistributedView(
        "myvar", vals, 1, { "time", "n" }, FileNetCDF::TimeRange{ 0, (std::size_t)ntimes });
    }
    {
      FileNetCDF in(comm, "test_netcdf_view_packed_time_chunk.nc", FILE_READONLY);
      Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::HostSpace> vals_in(
        "vals_in", ntimes, nlocal);
      in.readDistributedView(
        "myvar", vals_in, 1, { "time", "n" }, FileNetCDF::TimeRange{ 0, (std::size_t)ntimes });
      for (int t = 0; t != ntimes; ++t)
        for (int i = 0; i != nlocal; ++i) CHECK_CLOSE(vals(t, i), vals_in(t, i), 1.e-10);
    }
  }

  TEST(VIEW_WRITEREAD_PACKED_TIME_INDEX_1DOF)
  {
    // shape [N] at a single time_index out of a [NTIMES, N] dataset:
    // Kokkos::View<T*>, decomposed_axis=0 (one less than the chunked form's 1,
    // since the leading NTIMES axis is dropped).
    auto comm = getDefaultComm();
    int rank = comm->getRank();
    int nlocal = 3;
    int ntimes = 4;
    int offset = rank * nlocal;

    {
      FileNetCDF out(comm, "test_netcdf_view_packed_time_index.nc", FILE_CREATE);
      for (int t = 0; t != ntimes; ++t) {
        Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> vals("vals", nlocal);
        for (int i = 0; i != nlocal; ++i) vals(i) = 1000 * t + 10 * (offset + i);
        out.writeDistributedView("myvar", vals, 0, { "time", "n" }, (std::size_t)t);
      }
    }
    {
      FileNetCDF in(comm, "test_netcdf_view_packed_time_index.nc", FILE_READONLY);
      for (int t = 0; t != ntimes; ++t) {
        Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> vals_in("vals_in", nlocal);
        in.readDistributedView("myvar", vals_in, 0, { "time", "n" }, (std::size_t)t);
        for (int i = 0; i != nlocal; ++i)
          CHECK_CLOSE(1000 * t + 10 * (offset + i), vals_in(i), 1.e-10);
      }
    }
  }

  TEST(VIEW_WRITEREAD_PACKED_TIME_CHUNK_MULTIDOF)
  {
    // shape [NTIMES, NDOFS, N]: Kokkos::View<T***>, decomposed_axis=2, with time_chunk.
    auto comm = getDefaultComm();
    int rank = comm->getRank();
    int nlocal = 3;
    int ndofs = 2;
    int ntimes = 4;
    int offset = rank * nlocal;

    Kokkos::View<double***, Kokkos::LayoutRight, Kokkos::HostSpace> vals(
      "vals", ntimes, ndofs, nlocal);
    for (int t = 0; t != ntimes; ++t)
      for (int j = 0; j != ndofs; ++j)
        for (int i = 0; i != nlocal; ++i) vals(t, j, i) = 1000 * t + 100 * (offset + i) + j;

    {
      FileNetCDF out(comm, "test_netcdf_view_packed_time_chunk_multidof.nc", FILE_CREATE);
      out.writeDistributedView("myvar",
                               vals,
                               2,
                               { "time", "ndofs", "n" },
                               FileNetCDF::TimeRange{ 0, (std::size_t)ntimes });
    }
    {
      FileNetCDF in(comm, "test_netcdf_view_packed_time_chunk_multidof.nc", FILE_READONLY);
      Kokkos::View<double***, Kokkos::LayoutRight, Kokkos::HostSpace> vals_in(
        "vals_in", ntimes, ndofs, nlocal);
      in.readDistributedView("myvar",
                             vals_in,
                             2,
                             { "time", "ndofs", "n" },
                             FileNetCDF::TimeRange{ 0, (std::size_t)ntimes });
      for (int t = 0; t != ntimes; ++t)
        for (int j = 0; j != ndofs; ++j)
          for (int i = 0; i != nlocal; ++i) CHECK_CLOSE(vals(t, j, i), vals_in(t, j, i), 1.e-10);
    }
  }

  TEST(VIEW_WRITEREAD_PACKED_TIME_INDEX_MULTIDOF)
  {
    // shape [NDOFS, N] at a single time_index out of a [NTIMES, NDOFS, N]
    // dataset: Kokkos::View<T**>, decomposed_axis=1 (one less than the chunked
    // form's 2).
    auto comm = getDefaultComm();
    int rank = comm->getRank();
    int nlocal = 3;
    int ndofs = 2;
    int ntimes = 3;
    int offset = rank * nlocal;

    {
      FileNetCDF out(comm, "test_netcdf_view_packed_time_index_multidof.nc", FILE_CREATE);
      for (int t = 0; t != ntimes; ++t) {
        Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::HostSpace> vals("vals", ndofs, nlocal);
        for (int j = 0; j != ndofs; ++j)
          for (int i = 0; i != nlocal; ++i) vals(j, i) = 1000 * t + 100 * (offset + i) + j;
        out.writeDistributedView("myvar", vals, 1, { "time", "ndofs", "n" }, (std::size_t)t);
      }
    }
    {
      FileNetCDF in(comm, "test_netcdf_view_packed_time_index_multidof.nc", FILE_READONLY);
      for (int t = 0; t != ntimes; ++t) {
        Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::HostSpace> vals_in(
          "vals_in", ndofs, nlocal);
        in.readDistributedView("myvar", vals_in, 1, { "time", "ndofs", "n" }, (std::size_t)t);
        for (int j = 0; j != ndofs; ++j)
          for (int i = 0; i != nlocal; ++i)
            CHECK_CLOSE(1000 * t + 100 * (offset + i) + j, vals_in(j, i), 1.e-10);
      }
    }
  }

  TEST(VIEW_WRITEREAD_INT)
  {
    // int scalar_type instantiation, shape [N], decomposed_axis=0.
    auto comm = getDefaultComm();
    int rank = comm->getRank();
    int nlocal = 4;
    int offset = rank * nlocal;

    Kokkos::View<int*, Kokkos::LayoutRight, Kokkos::HostSpace> vals("vals", nlocal);
    for (int i = 0; i != nlocal; ++i) vals(i) = 100 * (offset + i);

    {
      FileNetCDF out(comm, "test_netcdf_view_int.nc", FILE_CREATE);
      out.writeDistributedView("myvar", vals, 0, { "n" });
    }
    {
      FileNetCDF in(comm, "test_netcdf_view_int.nc", FILE_READONLY);
      Kokkos::View<int*, Kokkos::LayoutRight, Kokkos::HostSpace> vals_in("vals_in", nlocal);
      in.readDistributedView("myvar", vals_in, 0, { "n" });
      for (int i = 0; i != nlocal; ++i) CHECK_EQUAL(vals(i), vals_in(i));
    }
  }

  //
  // Columnar shapes: NCOLUMNS and NCELLS_PER_COLUMN are kept as two real,
  // separate axes (NCELLS_PER_COLUMN always innermost/fastest-varying,
  // local in full on every rank); decomposed_axis names NCOLUMNS's
  // position, never the last axis.
  //
  TEST(VIEW_WRITEREAD_COLUMNAR)
  {
    // shape [NCOLUMNS, NCELLS_PER_COLUMN]: Kokkos::View<T**>, decomposed_axis=0.
    auto comm = getDefaultComm();
    int rank = comm->getRank();
    int ncolumns_local = 4;
    int ncells_per_column = 3;
    int offset = rank * ncolumns_local;

    Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::HostSpace> vals(
      "vals", ncolumns_local, ncells_per_column);
    for (int c = 0; c != ncolumns_local; ++c)
      for (int k = 0; k != ncells_per_column; ++k) vals(c, k) = 100 * (offset + c) + k;

    {
      FileNetCDF out(comm, "test_netcdf_view_columnar.nc", FILE_CREATE);
      out.writeDistributedView("myvar", vals, 0, { "ncolumn", "ncell_per_column" });
    }
    {
      FileNetCDF in(comm, "test_netcdf_view_columnar.nc", FILE_READONLY);
      Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::HostSpace> vals_in(
        "vals_in", ncolumns_local, ncells_per_column);
      in.readDistributedView("myvar", vals_in, 0, { "ncolumn", "ncell_per_column" });
      for (int c = 0; c != ncolumns_local; ++c)
        for (int k = 0; k != ncells_per_column; ++k) CHECK_CLOSE(vals(c, k), vals_in(c, k), 1.e-10);
    }
  }

  TEST(VIEW_WRITEREAD_COLUMNAR_MULTIDOF)
  {
    // shape [NDOFS, NCOLUMNS, NCELLS_PER_COLUMN]: Kokkos::View<T***>, decomposed_axis=1.
    auto comm = getDefaultComm();
    int rank = comm->getRank();
    int ncolumns_local = 3;
    int ncells_per_column = 2;
    int ndofs = 2;
    int offset = rank * ncolumns_local;

    Kokkos::View<double***, Kokkos::LayoutRight, Kokkos::HostSpace> vals(
      "vals", ndofs, ncolumns_local, ncells_per_column);
    for (int j = 0; j != ndofs; ++j)
      for (int c = 0; c != ncolumns_local; ++c)
        for (int k = 0; k != ncells_per_column; ++k)
          vals(j, c, k) = 1000 * j + 100 * (offset + c) + k;

    {
      FileNetCDF out(comm, "test_netcdf_view_columnar_multidof.nc", FILE_CREATE);
      out.writeDistributedView("myvar", vals, 1, { "ndofs", "ncolumn", "ncell_per_column" });
    }
    {
      FileNetCDF in(comm, "test_netcdf_view_columnar_multidof.nc", FILE_READONLY);
      Kokkos::View<double***, Kokkos::LayoutRight, Kokkos::HostSpace> vals_in(
        "vals_in", ndofs, ncolumns_local, ncells_per_column);
      in.readDistributedView("myvar", vals_in, 1, { "ndofs", "ncolumn", "ncell_per_column" });
      for (int j = 0; j != ndofs; ++j)
        for (int c = 0; c != ncolumns_local; ++c)
          for (int k = 0; k != ncells_per_column; ++k)
            CHECK_CLOSE(vals(j, c, k), vals_in(j, c, k), 1.e-10);
    }
  }

  TEST(VIEW_WRITEREAD_COLUMNAR_TIME_CHUNK)
  {
    // shape [NTIMES, NCOLUMNS, NCELLS_PER_COLUMN]: Kokkos::View<T***>, decomposed_axis=1,
    // with time_chunk.
    auto comm = getDefaultComm();
    int rank = comm->getRank();
    int ncolumns_local = 3;
    int ncells_per_column = 2;
    int ntimes = 4;
    int offset = rank * ncolumns_local;

    Kokkos::View<double***, Kokkos::LayoutRight, Kokkos::HostSpace> vals(
      "vals", ntimes, ncolumns_local, ncells_per_column);
    for (int t = 0; t != ntimes; ++t)
      for (int c = 0; c != ncolumns_local; ++c)
        for (int k = 0; k != ncells_per_column; ++k)
          vals(t, c, k) = 1000 * t + 100 * (offset + c) + k;

    {
      FileNetCDF out(comm, "test_netcdf_view_columnar_time_chunk.nc", FILE_CREATE);
      out.writeDistributedView("myvar",
                               vals,
                               1,
                               { "time", "ncolumn", "ncell_per_column" },
                               FileNetCDF::TimeRange{ 0, (std::size_t)ntimes });
    }
    {
      FileNetCDF in(comm, "test_netcdf_view_columnar_time_chunk.nc", FILE_READONLY);
      Kokkos::View<double***, Kokkos::LayoutRight, Kokkos::HostSpace> vals_in(
        "vals_in", ntimes, ncolumns_local, ncells_per_column);
      in.readDistributedView("myvar",
                             vals_in,
                             1,
                             { "time", "ncolumn", "ncell_per_column" },
                             FileNetCDF::TimeRange{ 0, (std::size_t)ntimes });
      for (int t = 0; t != ntimes; ++t)
        for (int c = 0; c != ncolumns_local; ++c)
          for (int k = 0; k != ncells_per_column; ++k)
            CHECK_CLOSE(vals(t, c, k), vals_in(t, c, k), 1.e-10);
    }
  }

  TEST(VIEW_WRITEREAD_COLUMNAR_TIME_INDEX)
  {
    // shape [NCOLUMNS, NCELLS_PER_COLUMN] at a single time_index out of a
    // [NTIMES, NCOLUMNS, NCELLS_PER_COLUMN] dataset: Kokkos::View<T**>,
    // decomposed_axis=0 (one less than the chunked form's 1).
    auto comm = getDefaultComm();
    int rank = comm->getRank();
    int ncolumns_local = 3;
    int ncells_per_column = 2;
    int ntimes = 3;
    int offset = rank * ncolumns_local;

    {
      FileNetCDF out(comm, "test_netcdf_view_columnar_time_index.nc", FILE_CREATE);
      for (int t = 0; t != ntimes; ++t) {
        Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::HostSpace> vals(
          "vals", ncolumns_local, ncells_per_column);
        for (int c = 0; c != ncolumns_local; ++c)
          for (int k = 0; k != ncells_per_column; ++k)
            vals(c, k) = 1000 * t + 100 * (offset + c) + k;
        out.writeDistributedView(
          "myvar", vals, 0, { "time", "ncolumn", "ncell_per_column" }, (std::size_t)t);
      }
    }
    {
      FileNetCDF in(comm, "test_netcdf_view_columnar_time_index.nc", FILE_READONLY);
      for (int t = 0; t != ntimes; ++t) {
        Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::HostSpace> vals_in(
          "vals_in", ncolumns_local, ncells_per_column);
        in.readDistributedView(
          "myvar", vals_in, 0, { "time", "ncolumn", "ncell_per_column" }, (std::size_t)t);
        for (int c = 0; c != ncolumns_local; ++c)
          for (int k = 0; k != ncells_per_column; ++k)
            CHECK_CLOSE(1000 * t + 100 * (offset + c) + k, vals_in(c, k), 1.e-10);
      }
    }
  }

  TEST(VIEW_WRITEREAD_COLUMNAR_TIME_CHUNK_MULTIDOF)
  {
    // shape [NTIMES, NDOFS, NCOLUMNS, NCELLS_PER_COLUMN]: Kokkos::View<T****>,
    // decomposed_axis=2, with time_chunk.
    auto comm = getDefaultComm();
    int rank = comm->getRank();
    int ncolumns_local = 2;
    int ncells_per_column = 2;
    int ndofs = 2;
    int ntimes = 3;
    int offset = rank * ncolumns_local;

    Kokkos::View<double****, Kokkos::LayoutRight, Kokkos::HostSpace> vals(
      "vals", ntimes, ndofs, ncolumns_local, ncells_per_column);
    for (int t = 0; t != ntimes; ++t)
      for (int j = 0; j != ndofs; ++j)
        for (int c = 0; c != ncolumns_local; ++c)
          for (int k = 0; k != ncells_per_column; ++k)
            vals(t, j, c, k) = 1000 * t + 100 * j + 10 * (offset + c) + k;

    {
      FileNetCDF out(comm, "test_netcdf_view_columnar_time_chunk_multidof.nc", FILE_CREATE);
      out.writeDistributedView("myvar",
                               vals,
                               2,
                               { "time", "ndofs", "ncolumn", "ncell_per_column" },
                               FileNetCDF::TimeRange{ 0, (std::size_t)ntimes });
    }
    {
      FileNetCDF in(comm, "test_netcdf_view_columnar_time_chunk_multidof.nc", FILE_READONLY);
      Kokkos::View<double****, Kokkos::LayoutRight, Kokkos::HostSpace> vals_in(
        "vals_in", ntimes, ndofs, ncolumns_local, ncells_per_column);
      in.readDistributedView("myvar",
                             vals_in,
                             2,
                             { "time", "ndofs", "ncolumn", "ncell_per_column" },
                             FileNetCDF::TimeRange{ 0, (std::size_t)ntimes });
      for (int t = 0; t != ntimes; ++t)
        for (int j = 0; j != ndofs; ++j)
          for (int c = 0; c != ncolumns_local; ++c)
            for (int k = 0; k != ncells_per_column; ++k)
              CHECK_CLOSE(vals(t, j, c, k), vals_in(t, j, c, k), 1.e-10);
    }
  }

  //
  // Open-ended appending: a dataset's time dimension must support being
  // extended indefinitely, in pieces, not just written once in full.
  //
  TEST(VIEW_MULTI_CHUNK_APPEND)
  {
    // Write the same [NTIMES, N] dataset in three separate writeDistributedView calls,
    // each extending past the previous call's end, then read the whole
    // thing back in one chunk.
    auto comm = getDefaultComm();
    int rank = comm->getRank();
    int nlocal = 3;
    int ntimes = 7;
    int offset = rank * nlocal;

    Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::HostSpace> vals("vals", ntimes, nlocal);
    for (int t = 0; t != ntimes; ++t)
      for (int i = 0; i != nlocal; ++i) vals(t, i) = 1000 * t + 10 * (offset + i);

    {
      FileNetCDF out(comm, "test_netcdf_multi_chunk_append.nc", FILE_CREATE);
      // three pieces: [0,2), [2,5), [5,7) -- uneven sizes on purpose.
      std::vector<FileNetCDF::TimeRange> pieces{ { 0, 2 }, { 2, 5 }, { 5, 7 } };
      for (const auto& piece : pieces) {
        int count = piece.second - piece.first;
        Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::HostSpace> slice(
          "slice", count, nlocal);
        for (int t = 0; t != count; ++t)
          for (int i = 0; i != nlocal; ++i) slice(t, i) = vals(piece.first + t, i);
        out.writeDistributedView("myvar", slice, 1, { "time", "n" }, piece);
      }
    }

    {
      FileNetCDF in(comm, "test_netcdf_multi_chunk_append.nc", FILE_READONLY);
      Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::HostSpace> vals_in(
        "vals_in", ntimes, nlocal);
      in.readDistributedView(
        "myvar", vals_in, 1, { "time", "n" }, FileNetCDF::TimeRange{ 0, (std::size_t)ntimes });
      for (int t = 0; t != ntimes; ++t)
        for (int i = 0; i != nlocal; ++i) CHECK_CLOSE(vals(t, i), vals_in(t, i), 1.e-10);
    }
  }

  TEST(VIEW_WRITE_CHUNK_READ_INDEX)
  {
    // Write as one chunk, then read back one time index at a time --
    // proves the two granularities interoperate on the same dataset.
    auto comm = getDefaultComm();
    int rank = comm->getRank();
    int nlocal = 3;
    int ntimes = 4;
    int offset = rank * nlocal;

    Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::HostSpace> vals("vals", ntimes, nlocal);
    for (int t = 0; t != ntimes; ++t)
      for (int i = 0; i != nlocal; ++i) vals(t, i) = 1000 * t + 10 * (offset + i);

    {
      FileNetCDF out(comm, "test_netcdf_write_chunk_read_index.nc", FILE_CREATE);
      out.writeDistributedView(
        "myvar", vals, 1, { "time", "n" }, FileNetCDF::TimeRange{ 0, (std::size_t)ntimes });
    }

    {
      FileNetCDF in(comm, "test_netcdf_write_chunk_read_index.nc", FILE_READONLY);
      for (int t = 0; t != ntimes; ++t) {
        Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> s_in("s_in", nlocal);
        in.readDistributedView("myvar", s_in, 0, { "time", "n" }, (std::size_t)t);
        for (int i = 0; i != nlocal; ++i) CHECK_CLOSE(vals(t, i), s_in(i), 1.e-10);
      }
    }
  }

  TEST(VIEW_WRITE_INDEX_READ_CHUNK)
  {
    // Write one time index at a time, then read the whole range back as
    // one chunk -- the reverse of the previous test.
    auto comm = getDefaultComm();
    int rank = comm->getRank();
    int nlocal = 3;
    int ntimes = 4;
    int offset = rank * nlocal;

    Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::HostSpace> vals("vals", ntimes, nlocal);
    for (int t = 0; t != ntimes; ++t)
      for (int i = 0; i != nlocal; ++i) vals(t, i) = 1000 * t + 10 * (offset + i);

    {
      FileNetCDF out(comm, "test_netcdf_write_index_read_chunk.nc", FILE_CREATE);
      for (int t = 0; t != ntimes; ++t) {
        Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> slice("slice", nlocal);
        for (int i = 0; i != nlocal; ++i) slice(i) = vals(t, i);
        out.writeDistributedView("myvar", slice, 0, { "time", "n" }, (std::size_t)t);
      }
    }

    {
      FileNetCDF in(comm, "test_netcdf_write_index_read_chunk.nc", FILE_READONLY);
      Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::HostSpace> vals_in(
        "vals_in", ntimes, nlocal);
      in.readDistributedView(
        "myvar", vals_in, 1, { "time", "n" }, FileNetCDF::TimeRange{ 0, (std::size_t)ntimes });
      for (int t = 0; t != ntimes; ++t)
        for (int i = 0; i != nlocal; ++i) CHECK_CLOSE(vals(t, i), vals_in(t, i), 1.e-10);
    }
  }

  TEST(VIEW_REOPEN_AND_APPEND)
  {
    // Simulates indefinite appending across separate runs: each "run" opens
    // a fresh FileNetCDF object in FILE_READWRITE mode (constructor opens,
    // destructor closes -- no object spans more than one append), writes
    // one more time index past the previous run's end, then closes.  A
    // final independent read confirms all appended indices survived every
    // open/close cycle.
    auto comm = getDefaultComm();
    int rank = comm->getRank();
    int nlocal = 3;
    int ntimes = 5;
    int offset = rank * nlocal;

    Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::HostSpace> vals("vals", ntimes, nlocal);
    for (int t = 0; t != ntimes; ++t)
      for (int i = 0; i != nlocal; ++i) vals(t, i) = 1000 * t + 10 * (offset + i);

    {
      FileNetCDF create(comm, "test_netcdf_reopen_append.nc", FILE_CREATE);
    }
    for (int t = 0; t != ntimes; ++t) {
      FileNetCDF out(comm, "test_netcdf_reopen_append.nc", FILE_READWRITE);
      Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> slice("slice", nlocal);
      for (int i = 0; i != nlocal; ++i) slice(i) = vals(t, i);
      out.writeDistributedView("myvar", slice, 0, { "time", "n" }, (std::size_t)t);
    }

    {
      FileNetCDF in(comm, "test_netcdf_reopen_append.nc", FILE_READONLY);
      for (int t = 0; t != ntimes; ++t) {
        Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> s_in("s_in", nlocal);
        in.readDistributedView("myvar", s_in, 0, { "time", "n" }, (std::size_t)t);
        for (int i = 0; i != nlocal; ++i) CHECK_CLOSE(vals(t, i), s_in(i), 1.e-10);
      }
    }
  }

  TEST(ATTR_TIME_MULTI_CHUNK_APPEND)
  {
    // Same open-ended-append concern, but for a time-dependent attribute
    // (e.g. "time" itself): write in two separate chunk calls, each
    // extending past the previous one, then read the whole thing back.
    int ntimes = 5;
    Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> times("times", ntimes);
    for (int t = 0; t != ntimes; ++t) times(t) = 10.0 * t;

    {
      FileNetCDF out(getDefaultComm(), "test_netcdf_attr_time_multi_chunk.nc", FILE_CREATE);
      Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> first("first", 2);
      for (int t = 0; t != 2; ++t) first(t) = times(t);
      out.writeAttribute("time", first, FileNetCDF::TimeRange{ 0, 2 });

      Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> rest("rest", ntimes - 2);
      for (int t = 0; t != ntimes - 2; ++t) rest(t) = times(2 + t);
      out.writeAttribute("time", rest, FileNetCDF::TimeRange{ 2, (std::size_t)ntimes });
    }
    {
      FileNetCDF in(getDefaultComm(), "test_netcdf_attr_time_multi_chunk.nc", FILE_READONLY);
      Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> times_in("times_in", ntimes);
      in.readAttribute("time", times_in, FileNetCDF::TimeRange{ 0, (std::size_t)ntimes });
      for (int t = 0; t != ntimes; ++t) CHECK_CLOSE(times(t), times_in(t), 1.e-10);
    }
  }

  TEST(ATTR_TIME_REOPEN_AND_APPEND)
  {
    // Same reopen-across-runs concern as VIEW_REOPEN_AND_APPEND, but for a
    // time-dependent attribute, appended one index at a time.
    int ntimes = 4;

    {
      FileNetCDF create(getDefaultComm(), "test_netcdf_attr_time_reopen_append.nc", FILE_CREATE);
    }
    for (int t = 0; t != ntimes; ++t) {
      FileNetCDF out(getDefaultComm(), "test_netcdf_attr_time_reopen_append.nc", FILE_READWRITE);
      out.writeAttribute("time", 10.0 * t, (std::size_t)t);
    }

    {
      FileNetCDF in(getDefaultComm(), "test_netcdf_attr_time_reopen_append.nc", FILE_READONLY);
      for (int t = 0; t != ntimes; ++t)
        CHECK_CLOSE(10.0 * t, in.readAttribute<double>("time", (std::size_t)t), 1.e-10);
    }
  }


  //
  // LocalView -- scalar, 1D, or 2D, identical on every rank, NOT
  // decomposed.  Mirrors FileHDF5's LOCALVIEW_* tests.
  //
  TEST(LOCALVIEW_WRITEREAD_SCALAR)
  {
    auto comm = getDefaultComm();
    double dval(3.14);
    int ival(42);

    {
      FileNetCDF out(comm, "test_netcdf_localview_scalar.nc", FILE_CREATE);
      out.writeLocalView("dval", dval);
      out.writeLocalView("ival", ival);
    }
    {
      FileNetCDF in(comm, "test_netcdf_localview_scalar.nc", FILE_READONLY);
      CHECK_CLOSE(dval, in.readLocalView<double>("dval"), 1.e-10);
      CHECK_EQUAL(ival, in.readLocalView<int>("ival"));
    }
  }

  TEST(LOCALVIEW_WRITEREAD_1D)
  {
    auto comm = getDefaultComm();
    int n = 6;
    Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> vals("vals", n);
    for (int i = 0; i != n; ++i) vals(i) = 10.0 * i;

    {
      FileNetCDF out(comm, "test_netcdf_localview_1d.nc", FILE_CREATE);
      out.writeLocalView("myvar", vals);
    }
    {
      FileNetCDF in(comm, "test_netcdf_localview_1d.nc", FILE_READONLY);
      Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> vals_in("vals_in", n);
      in.readLocalView("myvar", vals_in);
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
      FileNetCDF out(comm, "test_netcdf_localview_2d.nc", FILE_CREATE);
      out.writeLocalView("myvar", vals);
    }
    {
      FileNetCDF in(comm, "test_netcdf_localview_2d.nc", FILE_READONLY);
      Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::HostSpace> vals_in("vals_in", n, ndofs);
      in.readLocalView("myvar", vals_in);
      for (int i = 0; i != n; ++i)
        for (int j = 0; j != ndofs; ++j) CHECK_CLOSE(vals(i, j), vals_in(i, j), 1.e-10);
    }
  }

  // A default-constructed (0-size) view passed to readLocalView must be
  // resized to the on-disk shape rather than throwing a shape-mismatch
  // error -- mirrors FileHDF5's identical contract.
  TEST(LOCALVIEW_READ_1D_RESIZES_DEFAULT_CONSTRUCTED_VIEW)
  {
    auto comm = getDefaultComm();
    int n = 6;
    Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> vals("vals", n);
    for (int i = 0; i != n; ++i) vals(i) = 10.0 * i;

    {
      FileNetCDF out(comm, "test_netcdf_localview_1d_resize.nc", FILE_CREATE);
      out.writeLocalView("myvar", vals);
    }
    {
      FileNetCDF in(comm, "test_netcdf_localview_1d_resize.nc", FILE_READONLY);
      Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> vals_in; // default-constructed
      CHECK_EQUAL(0, (int)vals_in.extent(0));
      in.readLocalView("myvar", vals_in);
      CHECK_EQUAL(n, (int)vals_in.extent(0));
      for (int i = 0; i != n; ++i) CHECK_CLOSE(vals(i), vals_in(i), 1.e-10);
    }
  }

  // LocalView, paged by a single time_index -- on-disk array has one more
  // leading dimension than value.
  TEST(LOCALVIEW_WRITEREAD_TIME_INDEX_SCALAR)
  {
    auto comm = getDefaultComm();
    int ntimes = 4;

    {
      FileNetCDF out(comm, "test_netcdf_localview_time_index_scalar.nc", FILE_CREATE);
      for (std::size_t t = 0; t != (std::size_t)ntimes; ++t)
        out.writeLocalView("myvar", 100.0 * t, t);
    }
    {
      FileNetCDF in(comm, "test_netcdf_localview_time_index_scalar.nc", FILE_READONLY);
      for (std::size_t t = 0; t != (std::size_t)ntimes; ++t)
        CHECK_CLOSE(100.0 * t, in.readLocalView<double>("myvar", t), 1.e-10);
    }
  }

  TEST(LOCALVIEW_WRITEREAD_TIME_INDEX_1D)
  {
    auto comm = getDefaultComm();
    int n = 3, ntimes = 4;

    {
      FileNetCDF out(comm, "test_netcdf_localview_time_index_1d.nc", FILE_CREATE);
      for (int t = 0; t != ntimes; ++t) {
        Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> vals("vals", n);
        for (int i = 0; i != n; ++i) vals(i) = 1000 * t + i;
        out.writeLocalView("myvar", vals, (std::size_t)t);
      }
    }
    {
      FileNetCDF in(comm, "test_netcdf_localview_time_index_1d.nc", FILE_READONLY);
      for (int t = 0; t != ntimes; ++t) {
        Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> vals_in("vals_in", n);
        in.readLocalView("myvar", vals_in, (std::size_t)t);
        for (int i = 0; i != n; ++i) CHECK_CLOSE(1000 * t + i, vals_in(i), 1.e-10);
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
      FileNetCDF out(comm, "test_netcdf_localview_time_chunk_1d.nc", FILE_CREATE);
      out.writeLocalView("myvar", buf, FileNetCDF::TimeRange{ 0, (std::size_t)ntimes });
    }
    {
      FileNetCDF in(comm, "test_netcdf_localview_time_chunk_1d.nc", FILE_READONLY);
      Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::HostSpace> buf_in("buf_in", ntimes, n);
      in.readLocalView("myvar", buf_in, FileNetCDF::TimeRange{ 0, (std::size_t)ntimes });
      for (int t = 0; t != ntimes; ++t)
        for (int i = 0; i != n; ++i) CHECK_CLOSE(buf(t, i), buf_in(t, i), 1.e-10);
    }
  }

  // Open-ended appending: a LocalView dataset's time dimension must
  // support being extended indefinitely, in pieces.
  TEST(LOCALVIEW_MULTI_CHUNK_APPEND)
  {
    auto comm = getDefaultComm();
    int n = 3, ntimes = 7;

    Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::HostSpace> buf("buf", ntimes, n);
    for (int t = 0; t != ntimes; ++t)
      for (int i = 0; i != n; ++i) buf(t, i) = 1000 * t + i;

    {
      FileNetCDF out(comm, "test_netcdf_localview_multi_chunk_append.nc", FILE_CREATE);
      std::vector<FileNetCDF::TimeRange> pieces{ { 0, 2 }, { 2, 5 }, { 5, 7 } };
      for (const auto& piece : pieces) {
        int count = piece.second - piece.first;
        Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::HostSpace> slice("slice", count, n);
        for (int t = 0; t != count; ++t)
          for (int i = 0; i != n; ++i) slice(t, i) = buf(piece.first + t, i);
        out.writeLocalView("myvar", slice, piece);
      }
    }

    Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::HostSpace> buf_in("buf_in", ntimes, n);
    {
      FileNetCDF in(comm, "test_netcdf_localview_multi_chunk_append.nc", FILE_READONLY);
      in.readLocalView("myvar", buf_in, FileNetCDF::TimeRange{ 0, (std::size_t)ntimes });
    }
    for (int t = 0; t != ntimes; ++t)
      for (int i = 0; i != n; ++i) CHECK_CLOSE(buf(t, i), buf_in(t, i), 1.e-10);
  }

  // Error path -- a non-contiguous LayoutRight view must throw
  // DBC::Assertion rather than silently reading/writing a wrong buffer;
  // mirrors FileHDF5's identical check.
  TEST(WRITE_LOCALVIEW_2D_NONCONTIGUOUS_THROWS)
  {
    auto comm = getDefaultComm();
    int n = 3, ndofs = 2;
    Kokkos::LayoutRight padded_layout(n, ndofs);
    padded_layout.stride = ndofs + 2; // pad past the true row length
    Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::HostSpace> vals("vals", padded_layout);
    CHECK(!vals.span_is_contiguous());

    FileNetCDF out(comm, "test_netcdf_localview_2d_noncontiguous.nc", FILE_CREATE);
    CHECK_THROW(out.writeLocalView("myvar", vals), DBC::Assertion);
  }

  // Remaining error paths -- mirrors FileHDF5's *_THROWS coverage, adapted
  // to FileNetCDF's actual checks: readLocalView explicitly validates
  // on-disk rank/shape (see FileNetCDF_impl.hh), exactly like FileHDF5's
  // readLocalView; readDistributedView/readVector (the decomposed/DistributedView
  // analogue), by contrast, has no such explicit check at all -- it simply
  // hands the caller's view's own extents to nc_get_vara_*, which only
  // throws (via checkErr_) if the resulting start+count falls outside a
  // fixed dimension's bounds.  That's a real, intentional gap relative to
  // FileHDF5::readDistributedView (which does check), not an oversight
  // here -- so the readDistributedView tests below document the weaker guarantee
  // rather than assert a check that doesn't exist.
  TEST(READ_LOCALVIEW_1D_WRONG_LENGTH_THROWS)
  {
    auto comm = getDefaultComm();
    Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> vals("vals", 5);
    for (int i = 0; i != 5; ++i) vals(i) = i;

    {
      FileNetCDF out(comm, "test_netcdf_localview_1d_wrong_length.nc", FILE_CREATE);
      out.writeLocalView("myvar", vals);
    }
    {
      FileNetCDF in(comm, "test_netcdf_localview_1d_wrong_length.nc", FILE_READONLY);
      Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> vals_in(
        "vals_in", 4); // on-disk is length 5
      CHECK_THROW(in.readLocalView("myvar", vals_in), Errors::Message);
    }
  }

  TEST(READ_LOCALVIEW_RANK_MISMATCH_THROWS)
  {
    auto comm = getDefaultComm();
    Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> vals("vals", 6);
    for (int i = 0; i != 6; ++i) vals(i) = i;

    {
      FileNetCDF out(comm, "test_netcdf_localview_rank_mismatch.nc", FILE_CREATE);
      out.writeLocalView("myvar", vals); // on-disk rank 1
    }
    {
      FileNetCDF in(comm, "test_netcdf_localview_rank_mismatch.nc", FILE_READONLY);
      Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::HostSpace> vals_in(
        "vals_in", 2, 3); // same span, rank 2
      CHECK_THROW(in.readLocalView("myvar", vals_in), Errors::Message);
    }
  }

  TEST(READ_LOCALVIEW_2D_WRONG_SHAPE_THROWS)
  {
    auto comm = getDefaultComm();
    Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::HostSpace> vals("vals", 4, 3);
    for (int i = 0; i != 4; ++i)
      for (int j = 0; j != 3; ++j) vals(i, j) = 10 * i + j;

    {
      FileNetCDF out(comm, "test_netcdf_localview_2d_wrong_shape.nc", FILE_CREATE);
      out.writeLocalView("myvar", vals); // on-disk [4, 3]
    }
    {
      FileNetCDF in(comm, "test_netcdf_localview_2d_wrong_shape.nc", FILE_READONLY);
      Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::HostSpace> vals_in(
        "vals_in", 4, 2); // dim 1 mismatch
      CHECK_THROW(in.readLocalView("myvar", vals_in), Errors::Message);
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
      FileNetCDF out(comm, "test_netcdf_localview_2d_transposed.nc", FILE_CREATE);
      out.writeLocalView("myvar", vals); // on-disk [5, 2]
    }
    {
      FileNetCDF in(comm, "test_netcdf_localview_2d_transposed.nc", FILE_READONLY);
      Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::HostSpace> vals_in(
        "vals_in", ndofs, n); // [2, 5] -- transposed, same span
      CHECK_THROW(in.readLocalView("myvar", vals_in), Errors::Message);
    }
  }

  TEST(READ_LOCALVIEW_1D_TOO_FEW_ENTRIES_THROWS)
  {
    auto comm = getDefaultComm();
    Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> vals("vals", 8);
    for (int i = 0; i != 8; ++i) vals(i) = i;

    {
      FileNetCDF out(comm, "test_netcdf_localview_1d_too_few.nc", FILE_CREATE);
      out.writeLocalView("myvar", vals); // on-disk length 8
    }
    {
      FileNetCDF in(comm, "test_netcdf_localview_1d_too_few.nc", FILE_READONLY);
      Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> vals_in("vals_in", 3);
      CHECK_THROW(in.readLocalView("myvar", vals_in), Errors::Message);
    }
  }

  TEST(READ_LOCALVIEW_1D_TOO_MANY_ENTRIES_THROWS)
  {
    auto comm = getDefaultComm();
    Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> vals("vals", 8);
    for (int i = 0; i != 8; ++i) vals(i) = i;

    {
      FileNetCDF out(comm, "test_netcdf_localview_1d_too_many.nc", FILE_CREATE);
      out.writeLocalView("myvar", vals); // on-disk length 8
    }
    {
      FileNetCDF in(comm, "test_netcdf_localview_1d_too_many.nc", FILE_READONLY);
      Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> vals_in("vals_in", 20);
      CHECK_THROW(in.readLocalView("myvar", vals_in), Errors::Message);
    }
  }

  TEST(READ_LOCALVIEW_TIME_INDEX_WRONG_SHAPE_THROWS)
  {
    auto comm = getDefaultComm();
    Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> vals("vals", 3);
    for (int i = 0; i != 3; ++i) vals(i) = i;

    {
      FileNetCDF out(comm, "test_netcdf_localview_time_index_wrong_shape.nc", FILE_CREATE);
      out.writeLocalView("myvar", vals, (std::size_t)0); // on-disk [1, 3]
    }
    {
      FileNetCDF in(comm, "test_netcdf_localview_time_index_wrong_shape.nc", FILE_READONLY);
      Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> vals_in("vals_in",
                                                                            4); // dim 1 mismatch
      CHECK_THROW(in.readLocalView("myvar", vals_in, (std::size_t)0), Errors::Message);
    }
  }

  TEST(READ_LOCALVIEW_TIME_CHUNK_WRONG_SHAPE_THROWS)
  {
    auto comm = getDefaultComm();
    Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::HostSpace> buf("buf", 4, 3);
    for (int t = 0; t != 4; ++t)
      for (int i = 0; i != 3; ++i) buf(t, i) = 10 * t + i;

    {
      FileNetCDF out(comm, "test_netcdf_localview_time_chunk_wrong_shape.nc", FILE_CREATE);
      out.writeLocalView("myvar", buf, FileNetCDF::TimeRange{ 0, 4 }); // on-disk [4, 3]
    }
    {
      FileNetCDF in(comm, "test_netcdf_localview_time_chunk_wrong_shape.nc", FILE_READONLY);
      Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::HostSpace> buf_in(
        "buf_in", 4, 2); // dim 1 mismatch
      CHECK_THROW(in.readLocalView("myvar", buf_in, (FileNetCDF::TimeRange{ 0, 4 })),
                  Errors::Message);
    }
  }

  // readDistributedView/readVector (decomposed path) must throw when a rank's
  // requested slice runs past the on-disk "ncolumn" dimension's actual
  // (fixed, non-unlimited) length -- the equivalent of FileHDF5::
  // readDistributedView's explicit global-length check.  As of this
  // writing this does NOT throw; see the bug report sent to the
  // FileNetCDF implementer.  Left asserting the correct (currently
  // failing) behavior rather than weakened to match the bug.
  TEST(READ_VIEW_1D_PAST_GLOBAL_LENGTH_THROWS)
  {
    auto comm = getDefaultComm();
    int rank = comm->getRank();
    int nlocal = 5;
    Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> vals("vals", nlocal);
    for (int i = 0; i != nlocal; ++i) vals(i) = 10 * rank + i;

    {
      FileNetCDF out(comm, "test_netcdf_view_1d_past_length.nc", FILE_CREATE);
      out.writeDistributedView("myvar", vals, 0, { "n" }); // on-disk global length nlocal * nranks
    }
    {
      FileNetCDF in(comm, "test_netcdf_view_1d_past_length.nc", FILE_READONLY);
      // every rank asks for nlocal+1 entries at the same offset it wrote
      // from -- on the last rank this reads one entry past the fixed "n"
      // dimension's actual length.
      Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> vals_in("vals_in", nlocal + 1);
      CHECK_THROW(in.readDistributedView("myvar", vals_in, 0, { "n" }), Errors::Message);
    }
  }

  // Error paths -- non-contiguous views, for the LocalView overloads not
  // yet covered above (time_index-paged, time_chunk-paged, and the
  // decomposed writeDistributedView/readDistributedView path).
  TEST(WRITE_LOCALVIEW_TIME_INDEX_NONCONTIGUOUS_THROWS)
  {
    auto comm = getDefaultComm();
    int n = 3, ndofs = 2;
    Kokkos::LayoutRight padded_layout(n, ndofs);
    padded_layout.stride = ndofs + 2;
    Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::HostSpace> vals("vals", padded_layout);
    CHECK(!vals.span_is_contiguous());

    FileNetCDF out(comm, "test_netcdf_localview_time_index_noncontiguous.nc", FILE_CREATE);
    CHECK_THROW(out.writeLocalView("myvar", vals, (std::size_t)0), DBC::Assertion);
  }

  TEST(READ_LOCALVIEW_2D_NONCONTIGUOUS_THROWS)
  {
    auto comm = getDefaultComm();
    int n = 3, ndofs = 2;
    Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::HostSpace> vals("vals", n, ndofs);
    for (int i = 0; i != n; ++i)
      for (int j = 0; j != ndofs; ++j) vals(i, j) = 10 * i + j;

    {
      FileNetCDF out(comm, "test_netcdf_localview_2d_noncontiguous_read.nc", FILE_CREATE);
      out.writeLocalView("myvar", vals);
    }

    Kokkos::LayoutRight padded_layout(n, ndofs);
    padded_layout.stride = ndofs + 2;
    Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::HostSpace> vals_in("vals_in",
                                                                           padded_layout);
    CHECK(!vals_in.span_is_contiguous());

    FileNetCDF in(comm, "test_netcdf_localview_2d_noncontiguous_read.nc", FILE_READONLY);
    CHECK_THROW(in.readLocalView("myvar", vals_in), DBC::Assertion);
  }

  TEST(WRITE_LOCALVIEW_TIME_CHUNK_NONCONTIGUOUS_THROWS)
  {
    auto comm = getDefaultComm();
    int n = 3, ndofs = 2, ntimes = 2;
    Kokkos::LayoutRight padded_layout(ntimes, n, ndofs);
    padded_layout.stride = n * ndofs + 2;
    Kokkos::View<double***, Kokkos::LayoutRight, Kokkos::HostSpace> vals("vals", padded_layout);
    CHECK(!vals.span_is_contiguous());

    FileNetCDF out(comm, "test_netcdf_localview_time_chunk_noncontiguous.nc", FILE_CREATE);
    CHECK_THROW(out.writeLocalView("myvar", vals, FileNetCDF::TimeRange{ 0, (std::size_t)ntimes }),
                DBC::Assertion);
  }

  TEST(WRITE_VIEW_2D_NONCONTIGUOUS_THROWS)
  {
    auto comm = getDefaultComm();
    int nlocal = 4, ndofs = 3;
    Kokkos::LayoutRight padded_layout(nlocal, ndofs);
    padded_layout.stride = ndofs + 1;
    Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::HostSpace> vals("vals", padded_layout);
    CHECK(!vals.span_is_contiguous());

    FileNetCDF out(comm, "test_netcdf_view_2d_noncontiguous.nc", FILE_CREATE);
    CHECK_THROW(out.writeDistributedView("myvar", vals, 1, { "ndofs", "n" }), DBC::Assertion);
  }


  //
  // Vector -- delegates to the 1D writeDistributedView/readDistributedView, decomposed_axis=0.
  // Mirrors FileHDF5's VECTOR_WRITEREAD test.
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
      FileNetCDF out(comm, "test_netcdf_vector.nc", FILE_CREATE);
      out.writeVector("myvec", vec, { "n" });
    }

    Vector_type vec_in(map);
    {
      FileNetCDF in(comm, "test_netcdf_vector.nc", FILE_READONLY);
      in.readVector("myvec", vec_in, { "n" });
    }

    vec_in.update(-1, vec, 1);
    CHECK_CLOSE(0.0, vec_in.norm2(), 1.e-10);
  }


  //
  // writeDistributedView/readDistributedView, rank-2 LayoutLeft, decomposed_axis=0 -- exactly a
  // Tpetra::MultiVector's local view ([nlocal, ndofs], column-major),
  // read/written directly with a single nc_put_vara/nc_get_vara call (no
  // transpose, no per-DoF loop).
  //
  TEST(VIEW_WRITEREAD_MULTIVECTOR_LAYOUTLEFT)
  {
    auto comm = getDefaultComm();
    int size = comm->getSize();
    int rank = comm->getRank();
    int ndofs = 3;
    Map_ptr_type map = Teuchos::rcp(new Map_type(size * 4, 4, 0, comm));

    MultiVector_type vec(map, ndofs);
    {
      auto vecv = vec.getLocalViewHost(Tpetra::Access::ReadWrite);
      for (int i = 0; i != (int)vecv.extent(0); ++i)
        for (int j = 0; j != ndofs; ++j) vecv(i, j) = 100 * (i + rank * 4) + j;
    }

    {
      FileNetCDF out(comm, "test_netcdf_multivector_layoutleft.nc", FILE_CREATE);
      auto vecv = vec.getLocalViewHost(Tpetra::Access::ReadOnly);
      out.writeDistributedView("myvec", vecv, { "ndofs", "n" });
    }

    MultiVector_type vec_in(map, ndofs);
    {
      FileNetCDF in(comm, "test_netcdf_multivector_layoutleft.nc", FILE_READONLY);
      auto vecv_in = vec_in.getLocalViewHost(Tpetra::Access::ReadWrite);
      in.readDistributedView("myvec", vecv_in, { "ndofs", "n" });
    }

    auto vecv = vec.getLocalViewHost(Tpetra::Access::ReadOnly);
    auto vecv_in = vec_in.getLocalViewHost(Tpetra::Access::ReadOnly);
    for (int i = 0; i != (int)vecv.extent(0); ++i)
      for (int j = 0; j != ndofs; ++j) CHECK_CLOSE(vecv(i, j), vecv_in(i, j), 1.e-10);
  }
}
