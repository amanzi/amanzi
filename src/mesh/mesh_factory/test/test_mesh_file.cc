/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      William Perkins
*/

//! <MISSING_ONELINE_DOCSTRING>

#include <iostream>
#include <UnitTest++.h>

#include <AmanziComm.hh>

#include "dbc.hh"
#include "../FileFormat.hh"
#include "../MeshException.hh"

#define BOGUS_TEST_FILE "test/bogus.exo"
#define EXODUS_TEST_FILE "test/hex_3x3x3_ss.exo"
#define NEMESIS_TEST_FILE "test/hex_10x10x10_ss.par"
#define MOAB_TEST_FILE "test/hex_3x3x3_ss_4P.h5m"


SUITE(MeshFileType)
{
  TEST(ExodusII)
  {
    auto comm = Amanzi::getDefaultComm();

    // EXODUS_TEST_FILE is macro defined by cmake
    std::string fname(EXODUS_TEST_FILE);

    Amanzi::AmanziMesh::FileFormat f;
    try {
      f = Amanzi::AmanziMesh::fileFormatFromFilename(*comm, fname);
    } catch (const Amanzi::AmanziMesh::Message& e) {
      throw e;
    }

    CHECK(f == Amanzi::AmanziMesh::FileFormat::EXODUS_II);
  }

  TEST(Nemesis)
  {
    auto comm = Amanzi::getDefaultComm();

    // NEMESIS_TEST_FILE is macro defined by cmake
    std::string fname(NEMESIS_TEST_FILE);

    Amanzi::AmanziMesh::FileFormat f;
    if (comm->getSize() > 1 && comm->getSize() <= 4) {
      int ierr[1];
      ierr[0] = 0;
      try {
        f = Amanzi::AmanziMesh::fileFormatFromFilename(*comm, fname);
      } catch (const Amanzi::AmanziMesh::Message& e) {
        throw e;
      }

      CHECK(f == Amanzi::AmanziMesh::FileFormat::NEMESIS);
    } else {
      CHECK_THROW(f = Amanzi::AmanziMesh::fileFormatFromFilename(*comm, fname),
                  Amanzi::AmanziMesh::FileMessage);
    }
  }

  TEST(MOABHD5)
  {
    auto comm = Amanzi::getDefaultComm();

    // MOAB_TEST_FILE is macro defined by cmake
    std::string fname(MOAB_TEST_FILE);

    Amanzi::AmanziMesh::FileFormat f;
    try {
      f = Amanzi::AmanziMesh::fileFormatFromFilename(*comm, fname);
    } catch (const Amanzi::AmanziMesh::Message& e) {
      throw e;
    }

    CHECK(f == Amanzi::AmanziMesh::FileFormat::MOAB_HDF5);
  }

  TEST(PathFailure)
  {
    auto comm = Amanzi::getDefaultComm();

    std::string fname("/some/bogus/path.exo");

    Amanzi::AmanziMesh::FileFormat f;

    CHECK_THROW(Amanzi::AmanziMesh::fileFormatFromFilename(*comm, fname),
                Amanzi::AmanziMesh::FileMessage);
  }

  TEST(MagicNumberFailure)
  {
    auto comm = Amanzi::getDefaultComm();

    std::string fname(BOGUS_TEST_FILE);

    Amanzi::AmanziMesh::FileFormat f;

    CHECK_THROW(Amanzi::AmanziMesh::fileFormatFromFilename(*comm, fname),
                Amanzi::AmanziMesh::FileMessage);
  }
}
