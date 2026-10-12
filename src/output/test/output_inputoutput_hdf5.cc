/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

//! Tests InputHDF5/OutputHDF5 (checkpoint/restart via FileHDF5).
/*
  Exercises the Input/Output interface layer on top of FileHDF5's new
  readView/writeView/readAttribute/writeAttribute API: scalar attributes
  (double, int, string) and Vector_type/IntVector_type round-trips through
  the Input/Output interface, pulling Tpetra's host-resident local view out
  of the vector and handing FileHDF5 a plain Kokkos::View, exactly as
  InputOutputHDF5.cc itself does.

  These are the Vector_type-shaped round-trip tests formerly covered
  directly against FileHDF5 in output_file_hdf5.cc's VECTOR_WRITE_READ (now
  removed there, since Vector_type/MultiVector_type are no longer part of
  FileHDF5's own API) -- they belong at this layer instead, since going from
  a Tpetra vector to a Kokkos::View (and handling the Tpetra host view's
  LayoutLeft vs. FileHDF5's LayoutRight) is exactly what InputOutputHDF5
  does, not something FileHDF5 itself does anymore.

  Also covers a path-encoded attribute name (e.g. "42/time"), matching
  OutputXDMF's real usage of writeAttribute to scope per-cycle attributes
  into their own group without a separate path argument.
*/

#include "UnitTest++.h"

#include <string>

#include "AmanziComm.hh"
#include "AmanziTypes.hh"
#include "AmanziVector.hh"
#include "InputOutputHDF5.hh"
#include "FileHDF5.hh"

SUITE(INPUTOUTPUT_HDF5)
{
  using namespace Amanzi;

  TEST(ATTR_SCALAR_WRITEREAD)
  {
    auto comm = getDefaultComm();
    double dval = 6.1;
    int ival = 7;
    std::string sval("hello");

    Teuchos::ParameterList plist("checkpoint");
    plist.set<std::string>("file name base", "test_inputoutput_hdf5_attr");

    {
      OutputHDF5 out(plist, comm);
      out.createTimestep(0.0, 0);

      Teuchos::ParameterList attrs_d("dval");
      out.write(attrs_d, dval);
      Teuchos::ParameterList attrs_i("ival");
      out.write(attrs_i, ival);
      Teuchos::ParameterList attrs_s("sval");
      out.write(attrs_s, sval);

      out.finalizeTimestep();
    }

    std::string filename = plist.get<std::string>("file name base") + "00000.h5";
    {
      InputHDF5 in(comm, filename);

      double dval_in;
      Teuchos::ParameterList attrs_d("dval");
      in.read(attrs_d, dval_in);
      CHECK_CLOSE(dval, dval_in, 1.e-10);

      int ival_in;
      Teuchos::ParameterList attrs_i("ival");
      in.read(attrs_i, ival_in);
      CHECK_EQUAL(ival, ival_in);

      std::string sval_in;
      Teuchos::ParameterList attrs_s("sval");
      in.read(attrs_s, sval_in);
      CHECK_EQUAL(sval, sval_in);
    }
  }


  TEST(VECTOR_WRITE_READ)
  {
    auto comm = getDefaultComm();
    int size = comm->getSize();
    int rank = comm->getRank();
    Map_ptr_type map = Teuchos::rcp(new Map_type(size * 3, 3, 0, comm));

    Vector_ptr_type vec = Teuchos::rcp(new Vector_type(map));
    {
      auto vecv = vec->getLocalViewHost(Tpetra::Access::ReadWrite);
      for (int i = 0; i != (int)vecv.extent(0); ++i) { vecv(i, 0) = i + rank * 3; }
    }

    Teuchos::ParameterList plist("checkpoint");
    plist.set<std::string>("file name base", "test_inputoutput_hdf5_vector");

    {
      OutputHDF5 out(plist, comm);
      out.createTimestep(0.0, 0);
      Teuchos::ParameterList attrs("myvec");
      out.write(attrs, *vec);
      out.finalizeTimestep();
    }

    std::string filename = plist.get<std::string>("file name base") + "00000.h5";
    Vector_type vec_in(map);
    {
      InputHDF5 in(comm, filename);
      Teuchos::ParameterList attrs("myvec");
      in.read(attrs, vec_in);
    }

    vec_in.update(-1, *vec, 1);
    CHECK_CLOSE(0.0, vec_in.norm2(), 1.e-10);
  }


  TEST(INTVECTOR_WRITE_READ)
  {
    auto comm = getDefaultComm();
    int size = comm->getSize();
    int rank = comm->getRank();
    Map_ptr_type map = Teuchos::rcp(new Map_type(size * 4, 4, 0, comm));

    IntVector_type vec(map);
    {
      auto vecv = vec.getLocalViewHost(Tpetra::Access::ReadWrite);
      for (int i = 0; i != (int)vecv.extent(0); ++i) { vecv(i, 0) = i + rank * 4; }
    }

    Teuchos::ParameterList plist("checkpoint");
    plist.set<std::string>("file name base", "test_inputoutput_hdf5_intvector");

    {
      OutputHDF5 out(plist, comm);
      out.createTimestep(0.0, 0);
      Teuchos::ParameterList attrs("myintvec");
      out.write(attrs, vec);
      out.finalizeTimestep();
    }

    std::string filename = plist.get<std::string>("file name base") + "00000.h5";
    IntVector_type vec_in(map);
    {
      InputHDF5 in(comm, filename);
      Teuchos::ParameterList attrs("myintvec");
      in.read(attrs, vec_in);
    }

    auto vecv = vec.getLocalViewHost(Tpetra::Access::ReadOnly);
    auto vecv_in = vec_in.getLocalViewHost(Tpetra::Access::ReadOnly);
    for (int i = 0; i != (int)vecv.extent(0); ++i) CHECK_EQUAL(vecv(i, 0), vecv_in(i, 0));
  }


  // Matches OutputXDMF's real usage: writeAttribute(attrs.name(), val) where
  // attrs.name() itself is "<cycle>/<name>" -- verifies FileHDF5 correctly
  // splits this into group "<cycle>/" + attribute "<name>", round-tripping
  // through FileHDF5 directly (below the Input/Output layer, since
  // Input/Output's own read/write(attrs, val) uses attrs.name() verbatim as
  // the attribute name passed to FileHDF5, exactly reproducing what
  // OutputXDMF does).
  TEST(ATTR_PATH_ENCODED_NAME)
  {
    auto comm = getDefaultComm();
    double time_42 = 4.2;
    double time_43 = 4.3;

    {
      FileHDF5 out(comm, "test_inputoutput_hdf5_cycle_attr.h5", FILE_CREATE);
      out.writeAttribute<double>("time", "/42", time_42);
      out.writeAttribute<double>("time", "/43", time_43);
    }
    {
      FileHDF5 in(comm, "test_inputoutput_hdf5_cycle_attr.h5", FILE_READONLY);
      CHECK_CLOSE(time_42, in.readAttribute<double>("time", "/42"), 1.e-10);
      CHECK_CLOSE(time_43, in.readAttribute<double>("time", "/43"), 1.e-10);
    }
  }


  // Exercises InputHDF5::readRange directly: writes a small multi-DoF,
  // multi-time-slice field as per-DoF/per-time-index datasets (matching
  // readRange's "name.<dof>/<time_index>" convention, see
  // InputOutputHDF5.cc), then reads it back both in one full-range call and
  // in several smaller page-sized chunks, confirming values match either way.
  TEST(READRANGE_ROUNDTRIP)
  {
    auto comm = getDefaultComm();
    int size = comm->getSize();
    int rank = comm->getRank();
    int nlocal = 3;
    Map_ptr_type map = Teuchos::rcp(new Map_type(size * nlocal, nlocal, 0, comm));

    int ntimes = 4;
    int ndofs = 2;

    std::string filename = "test_inputoutput_hdf5_readrange.h5";
    {
      FileHDF5 out(comm, filename, FILE_CREATE);
      for (int t = 0; t != ntimes; ++t) {
        for (int j = 0; j != ndofs; ++j) {
          Vector_type vec(map);
          auto vecv = vec.getLocalViewHost(Tpetra::Access::ReadWrite);
          for (int i = 0; i != (int)vecv.extent(0); ++i) {
            vecv(i, 0) = 1000 * t + 10 * j + i + rank * nlocal;
          }
          std::string path = "myvar." + std::to_string(j) + "/";
          out.writeVector(std::to_string(t), path, vec);
        }
      }
    }

    auto expected = [&](int t, int j, int i) { return 1000 * t + 10 * j + i + rank * nlocal; };

    InputHDF5 in(comm, filename);
    Teuchos::ParameterList attrs("myvar");

    // one call covering the whole range
    {
      Kokkos::View<double***, Kokkos::LayoutRight, Kokkos::HostSpace> buf(
        "buf", ntimes, ndofs, nlocal);
      in.readRange(attrs, 0, ntimes, buf);
      for (int t = 0; t != ntimes; ++t)
        for (int j = 0; j != ndofs; ++j)
          for (int i = 0; i != nlocal; ++i)
            CHECK_CLOSE(expected(t, j, i), buf(t, j, i), 1.e-10);
    }

    // several smaller chunks
    {
      int chunk = 2;
      for (int t0 = 0; t0 < ntimes; t0 += chunk) {
        int count = std::min(chunk, ntimes - t0);
        Kokkos::View<double***, Kokkos::LayoutRight, Kokkos::HostSpace> buf(
          "buf", count, ndofs, nlocal);
        in.readRange(attrs, t0, count, buf);
        for (int t = 0; t != count; ++t)
          for (int j = 0; j != ndofs; ++j)
            for (int i = 0; i != nlocal; ++i)
              CHECK_CLOSE(expected(t0 + t, j, i), buf(t, j, i), 1.e-10);
      }
    }
  }
}
