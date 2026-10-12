/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

//! Tests InputOutputVerticallyStructuredNetCDF (InputVerticallyStructuredNetCDF/OutputVerticallyStructuredNetCDF).
/*
  Exercises InputOutputVerticallyStructuredNetCDF end-to-end for the
  subfieldnames-present and subfieldnames-absent cases, confirming variable
  naming matches OutputUtils::getNames's convention exactly; a test
  confirming ColumnarMeshInfo errors clearly when mesh->columns does not
  exist; and a test confirming multiple timesteps are written as time
  slices into one on-disk file rather than one file per cycle.

  The positive-path tests (anything that actually writes/reads columnar
  data) load domain.lnd.oakharbor_transect.exo via the MSTK framework with
  `"partitioner" = "from exodus file"` and `"build columns" = true` -- a
  real FROM_EXO-partitioned mesh (5 columns, 15 cells per column).  The
  file's own partition coloring is for 2 ranks, so these tests run on 1 or
  2 ranks (skipped otherwise), exercising ColumnarMeshInfo's uniform-height
  and memory-contiguity checks against real mesh structure, not merely
  assumed.  Skipped entirely if MSTK is not enabled.
*/

#include "UnitTest++.h"

#include "AmanziComm.hh"
#include "AmanziTypes.hh"
#include "GeometricModel.hh"
#include "MeshFactory.hh"
#include "CompositeVector.hh"
#include "InputOutputVerticallyStructuredNetCDF.hh"

SUITE(INPUTOUTPUT_VERTICALLY_STRUCTURED_NETCDF)
{
  using namespace Amanzi;

  // domain.lnd.oakharbor_transect.exo has 5 columns, 15 cells per column.
  const int ncells_per_column = 15;

  Teuchos::RCP<AmanziMesh::Mesh> createFromExoColumnarMesh(const Comm_ptr_type& comm)
  {
    auto plist = Teuchos::rcp(new Teuchos::ParameterList());
    plist->set<std::string>("partitioner", "from exodus file");
    plist->set<bool>("build columns", true);

    AmanziMesh::MeshFactory fac(comm, Teuchos::null, plist);
    fac.set_preference(AmanziMesh::Preference{ AmanziMesh::Framework::MSTK });
    return fac.create("test/domain.lnd.oakharbor_transect.exo");
  }


  TEST(REJECT_MESH_WITHOUT_COLUMNS)
  {
    // a mesh that never called buildColumns() -- mesh->columns is null --
    // must be rejected clearly by ColumnarMeshInfo, independent of the
    // FROM_EXO fixture above.
    auto comm = getDefaultComm();
    if (comm->getSize() > 1) return; // SIMPLE framework is serial-only

    AmanziMesh::MeshFactory fac(comm);
    fac.set_preference(AmanziMesh::Preference{ AmanziMesh::Framework::SIMPLE });
    auto mesh = fac.create(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2, 2, 3);

    Teuchos::ParameterList plist("domain");
    plist.set<std::string>("file name base", "test_inputoutput_vsnetcdf_reject");

    CHECK_THROW(OutputVerticallyStructuredNetCDF(plist, comm, mesh), Errors::Message);
  }


  TEST(SUBFIELDNAMES_PRESENT_SPLITS_PER_DOF)
  {
    auto comm = getDefaultComm();
    if (!AmanziMesh::framework_enabled(AmanziMesh::Framework::MSTK)) return;
    if (comm->getSize() != 1 && comm->getSize() != 2) return;

    auto mesh = createFromExoColumnarMesh(comm);

    int ndofs = 2;
    auto map = mesh->getMap(AmanziMesh::Entity_kind::CELL, false);
    int ncolumns_local = map->getLocalNumElements() / ncells_per_column;

    MultiVector_type vec(map, ndofs);
    {
      auto vecv = vec.getLocalViewHost(Tpetra::Access::ReadWrite);
      for (int i = 0; i != (int)vecv.extent(0); ++i) {
        int global_id = (int)map->getGlobalElement(i);
        vecv(i, 0) = global_id;
        vecv(i, 1) = 100 + global_id;
      }
    }

    Teuchos::ParameterList attrs("tcc");
    attrs.set<Teuchos::Array<std::string>>("subfieldnames", { "H+", "OH-" });
    attrs.set<AmanziMesh::Entity_kind>("location", AmanziMesh::Entity_kind::CELL);

    Teuchos::ParameterList plist("domain");
    plist.set<std::string>("file name base", "test_inputoutput_vsnetcdf_split");

    std::string filename;
    {
      OutputVerticallyStructuredNetCDF out(plist, comm, mesh);
      out.createTimestep(0.0, 0);
      out.write(attrs, vec);
      out.finalizeTimestep();
      filename = out.getFilename(0);
    }

    // confirm the on-disk variable names match OutputUtils::getNames's
    // convention exactly: "tcc.H+" and "tcc.OH-", not a single "tcc" with a
    // real DoF axis.  Each split variable has no DoF axis at all --
    // [ntimes, ncolumns, ncells_per_column], decomposed_axis=1 -- matching
    // OutputVerticallyStructuredNetCDF::writeRavelled_'s single-DoF shape
    // (paged by the one time index written above).  Each rank reads back
    // only its own ncolumns_local columns, matching what writeRavelled_
    // wrote for this rank.
    {
      FileNetCDF in(comm, filename, FILE_READONLY);

      Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::HostSpace> buf0(
        "buf0", ncolumns_local, ncells_per_column);
      in.readDistributedView(
        "tcc.H+", buf0, 0, { "time", "ncolumn", "ncell_per_column" }, (std::size_t)0);
      Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::HostSpace> buf1(
        "buf1", ncolumns_local, ncells_per_column);
      in.readDistributedView(
        "tcc.OH-", buf1, 0, { "time", "ncolumn", "ncell_per_column" }, (std::size_t)0);

      auto vecv = vec.getLocalViewHost(Tpetra::Access::ReadOnly);
      for (int i = 0; i != (int)vecv.extent(0); ++i) {
        int col = i / ncells_per_column;
        int layer = i % ncells_per_column;
        CHECK_CLOSE(vecv(i, 0), buf0(col, layer), 1.e-10);
        CHECK_CLOSE(vecv(i, 1), buf1(col, layer), 1.e-10);
      }
    }

    // InputVerticallyStructuredNetCDF only implements readRange() -- every
    // other Input method, including read(MultiVector_type&), is not
    // implemented and throws clearly.
    {
      InputVerticallyStructuredNetCDF in(comm, filename, mesh);
      MultiVector_type vec_in(map, ndofs);
      CHECK_THROW(in.read(attrs, vec_in), Errors::Message);
    }
  }


  // A surface mesh extracted (and flattened) from a 3D parent is still
  // embedded in 3-space -- getSpaceDimension() == 3 -- even though it is
  // topologically 2D (getManifoldDimension() == 2).  ColumnarMeshInfo must
  // key off manifold dimension, not space dimension, or it wrongly demands
  // mesh->columns on a mesh that can never have them.  A CELL field on
  // such a mesh must be written flat/unraveled (named "ncell"), not
  // columnar.  Mirrors ats_mesh_factory.cc's createMeshSurface, which
  // extracts the "surface" labeled set (label 2) from the same fixture
  // used by createFromExoColumnarMesh above.
  TEST(SURFACE_MESH_CELL_FIELD_IS_FLAT)
  {
    auto comm = getDefaultComm();
    if (!AmanziMesh::framework_enabled(AmanziMesh::Framework::MSTK)) return;
    if (comm->getSize() != 1 && comm->getSize() != 2) return;

    std::string exo_file = "test/domain.lnd.oakharbor_transect.exo";

    Teuchos::ParameterList reg_spec;
    auto& surf_reg_spec = reg_spec.sublist("surface").sublist("region: labeled set");
    surf_reg_spec.set<std::string>("label", "2");
    surf_reg_spec.set<std::string>("file", exo_file);
    surf_reg_spec.set<std::string>("format", "Exodus II");
    surf_reg_spec.set<std::string>("entity", "face");
    auto gm = Teuchos::rcp(new AmanziGeometry::GeometricModel(3, reg_spec, *comm));

    auto parent_plist = Teuchos::rcp(new Teuchos::ParameterList());
    parent_plist->set<std::string>("partitioner", "from exodus file");
    AmanziMesh::MeshFactory parent_fac(comm, gm, parent_plist);
    parent_fac.set_preference(AmanziMesh::Preference{ AmanziMesh::Framework::MSTK });
    auto parent_mesh = parent_fac.create(exo_file);
    CHECK_EQUAL(3, parent_mesh->getManifoldDimension());

    AmanziMesh::MeshFactory surf_fac(comm, gm, Teuchos::null);
    surf_fac.set_preference(AmanziMesh::Preference{ AmanziMesh::Framework::MSTK });
    std::vector<std::string> regions{ "surface" };
    // flatten=false: this submanifold mesh keeps its nodes' real 3D
    // coordinates (getSpaceDimension() stays 3) despite being topologically
    // 2D (getManifoldDimension() == 2) -- the actual bug.  (flatten=true
    // instead drops to a true 2D embedding, getSpaceDimension() == 2 too,
    // and would not exercise this.)
    auto mesh = surf_fac.create(parent_mesh, regions, AmanziMesh::Entity_kind::FACE, false);

    CHECK_EQUAL(3, mesh->getSpaceDimension());
    CHECK_EQUAL(2, mesh->getManifoldDimension());

    auto map = mesh->getMap(AmanziMesh::Entity_kind::CELL, false);
    Vector_type vec(map);
    {
      auto vecv = vec.getLocalViewHost(Tpetra::Access::ReadWrite);
      for (int i = 0; i != (int)vecv.extent(0); ++i) vecv(i, 0) = (double)map->getGlobalElement(i);
    }

    Teuchos::ParameterList attrs("surface-pressure");
    attrs.set<AmanziMesh::Entity_kind>("location", AmanziMesh::Entity_kind::CELL);

    Teuchos::ParameterList plist("domain");
    plist.set<std::string>("file name base", "test_inputoutput_vsnetcdf_surface");

    std::string filename;
    {
      OutputVerticallyStructuredNetCDF out(plist, comm, mesh);
      out.createTimestep(0.0, 0);
      out.write(attrs, vec);
      out.finalizeTimestep();
      filename = out.getFilename(0);
    }

    // confirm the on-disk variable is flat ([ncell], decomposed_axis=0),
    // not columnar -- no "ncolumn"/"ncell_per_column" dimension exists.
    {
      FileNetCDF in(comm, filename, FILE_READONLY);
      int n = (int)map->getLocalNumElements();
      Kokkos::View<double*, Kokkos::LayoutRight, Kokkos::HostSpace> buf("buf", n);
      in.readDistributedView(attrs.name(), buf, 0, { "time", "ncell" }, (std::size_t)0);

      auto vecv = vec.getLocalViewHost(Tpetra::Access::ReadOnly);
      for (int i = 0; i != n; ++i) CHECK_CLOSE(vecv(i, 0), buf(i), 1.e-10);
    }
  }


  // Confirms the actual bug fix this class exists for: writing several
  // timesteps through OutputVerticallyStructuredNetCDF produces ONE file
  // (getFilename does not vary with cycle) holding each cycle's data as a
  // distinct time slice, not one file per cycle -- and that each slice's
  // raw on-disk values match what was written.
  TEST(MULTIPLE_TIMESTEPS_SHARE_ONE_FILE)
  {
    auto comm = getDefaultComm();
    if (!AmanziMesh::framework_enabled(AmanziMesh::Framework::MSTK)) return;
    if (comm->getSize() != 1 && comm->getSize() != 2) return;

    auto mesh = createFromExoColumnarMesh(comm);
    auto map = mesh->getMap(AmanziMesh::Entity_kind::CELL, false);
    int ncolumns_local = map->getLocalNumElements() / ncells_per_column;
    int column_offset = map->getMinGlobalIndex() / ncells_per_column;
    Teuchos::ParameterList attrs("pressure");
    attrs.set<AmanziMesh::Entity_kind>("location", AmanziMesh::Entity_kind::CELL);

    Teuchos::ParameterList plist("domain");
    plist.set<std::string>("file name base", "test_inputoutput_vsnetcdf_multistep");

    int ntimes = 3;
    std::vector<std::string> filenames(ntimes);
    {
      OutputVerticallyStructuredNetCDF out(plist, comm, mesh);
      for (int t = 0; t != ntimes; ++t) {
        Vector_type vec(map);
        {
          auto vecv = vec.getLocalViewHost(Tpetra::Access::OverwriteAll);
          for (int i = 0; i != (int)vecv.extent(0); ++i)
            vecv(i, 0) = 1000 * t + map->getGlobalElement(i);
        }
        out.createTimestep((double)t, t);
        out.write(attrs, vec);
        out.finalizeTimestep();
        filenames[t] = out.getFilename(t);
      }
    }

    // every cycle reports the same, single filename
    for (int t = 1; t != ntimes; ++t) CHECK_EQUAL(filenames[0], filenames[t]);

    // and that one file holds all ntimes slices, each matching what was
    // written -- each rank reads back only its own ncolumns_local columns.
    FileNetCDF in(comm, filenames[0], FILE_READONLY);
    for (int t = 0; t != ntimes; ++t) {
      Kokkos::View<double**, Kokkos::LayoutRight, Kokkos::HostSpace> buf(
        "buf", ncolumns_local, ncells_per_column);
      in.readDistributedView(
        "pressure", buf, 0, { "time", "ncolumn", "ncell_per_column" }, (std::size_t)t);
      for (int col = 0; col != ncolumns_local; ++col) {
        for (int layer = 0; layer != ncells_per_column; ++layer) {
          int global_id = (column_offset + col) * ncells_per_column + layer;
          CHECK_CLOSE((double)(1000 * t + global_id), buf(col, layer), 1.e-10);
        }
      }
    }
  }


  // Exercises InputVerticallyStructuredNetCDF::readRange directly: writes a
  // small multi-DoF, multi-time-slice columnar field via raw FileNetCDF
  // calls (one writeDistributedView per time index, decomposed_axis=2,
  // matching the [NTIMES, NDOFS, NCOLUMNS, NCELLS_PER_COLUMN] convention
  // FileNetCDF.hh documents), then reads it back through
  // InputVerticallyStructuredNetCDF::readRange, both in one full-range call
  // and in several smaller page-sized chunks, confirming values match
  // either way and that the column/layer reshape (local_id = col *
  // ncells_per_column + layer) round-trips correctly.
  TEST(READRANGE_ROUNDTRIP)
  {
    auto comm = getDefaultComm();
    if (!AmanziMesh::framework_enabled(AmanziMesh::Framework::MSTK)) return;
    if (comm->getSize() != 1 && comm->getSize() != 2) return;

    auto mesh = createFromExoColumnarMesh(comm);
    auto map = mesh->getMap(AmanziMesh::Entity_kind::CELL, false);
    int ncolumns_local = map->getLocalNumElements() / ncells_per_column;
    int column_offset = map->getMinGlobalIndex() / ncells_per_column;

    int ndofs = 2;
    int ntimes = 4;
    std::string filename = "test_inputoutput_vsnetcdf_readrange.nc";

    // global_col/layer address the full [ncolumns, ncells_per_column] mesh;
    // expected is keyed by those global coordinates so every rank agrees
    // on what belongs where, independent of how many ranks are running.
    auto expected = [&](int t, int j, int global_col, int layer) {
      return 1000 * t + 100 * j + global_col * ncells_per_column + layer;
    };

    {
      FileNetCDF out(comm, filename, FILE_CREATE);
      for (int t = 0; t != ntimes; ++t) {
        // each rank writes only its own ncolumns_local columns -- writeDistributedView's
        // decomposed_axis (2, the ncolumns axis) does the global offset/extent
        // accounting via an MPI scan/reduceAll internally.
        Kokkos::View<double****, Kokkos::LayoutRight, Kokkos::HostSpace> buf(
          "buf", 1, ndofs, ncolumns_local, ncells_per_column);
        for (int j = 0; j != ndofs; ++j)
          for (int col = 0; col != ncolumns_local; ++col)
            for (int layer = 0; layer != ncells_per_column; ++layer)
              buf(0, j, col, layer) = expected(t, j, column_offset + col, layer);
        // dim_names must match what InputVerticallyStructuredNetCDF::readRange
        // itself requests below ("ndofs_velocity", not a generic "ndofs").
        out.writeDistributedView("velocity",
                                 buf,
                                 2,
                                 { "time", "ndofs_velocity", "ncolumn", "ncell_per_column" },
                                 FileNetCDF::TimeRange{ (std::size_t)t, (std::size_t)t + 1 });
      }
    }

    InputVerticallyStructuredNetCDF in(comm, filename, mesh);
    Teuchos::ParameterList attrs("velocity");
    int num_entities = ncolumns_local * ncells_per_column;

    // one call covering the whole range
    {
      Kokkos::View<double***, Kokkos::LayoutRight, Kokkos::HostSpace> buf(
        "buf", ntimes, ndofs, num_entities);
      in.readRange(attrs, 0, ntimes, buf);
      for (int t = 0; t != ntimes; ++t)
        for (int j = 0; j != ndofs; ++j)
          for (int col = 0; col != ncolumns_local; ++col)
            for (int layer = 0; layer != ncells_per_column; ++layer) {
              int local_id = col * ncells_per_column + layer;
              CHECK_CLOSE(expected(t, j, column_offset + col, layer), buf(t, j, local_id), 1.e-10);
            }
    }

    // several smaller chunks
    {
      int chunk = 2;
      for (int t0 = 0; t0 < ntimes; t0 += chunk) {
        int count = std::min(chunk, ntimes - t0);
        Kokkos::View<double***, Kokkos::LayoutRight, Kokkos::HostSpace> buf(
          "buf", count, ndofs, num_entities);
        in.readRange(attrs, t0, count, buf);
        for (int t = 0; t != count; ++t)
          for (int j = 0; j != ndofs; ++j)
            for (int col = 0; col != ncolumns_local; ++col)
              for (int layer = 0; layer != ncells_per_column; ++layer) {
                int local_id = col * ncells_per_column + layer;
                CHECK_CLOSE(
                  expected(t0 + t, j, column_offset + col, layer), buf(t, j, local_id), 1.e-10);
              }
      }
    }
  }
}
