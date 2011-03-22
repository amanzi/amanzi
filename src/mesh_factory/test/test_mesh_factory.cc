// -------------------------------------------------------------
// file: test_mesh_factory.cc
// -------------------------------------------------------------
// -------------------------------------------------------------
// Battelle Memorial Institute
// Pacific Northwest Laboratory
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created March 18, 2011 by William A. Perkins
// Last Change: Tue Mar 22 12:00:06 2011 by William A. Perkins <d3g096@PE10900.pnl.gov>
// -------------------------------------------------------------


static const char* SCCS_ID = "$Id$ Battelle PNL";

#include <iostream>
#include <UnitTest++.h>

#include <Epetra_MpiComm.h>

#include "dbc.hh"
#include "../MeshFactory.hh"
#include "../FrameworkTraits.hh"

// Check to see if we have some files to read

#ifndef BOGUS_TEST_FILE
#define BOGUS_TEST_FILE "bogus.file"
#endif

#ifndef EXODUS_TEST_FILE
#error EXODUS_TEST_FILE must be defined
#endif

#ifndef NEMESIS_TEST_FILE
#error EXODUS_TEST_FILE must be defined
#endif

#ifndef MOAB_TEST_FILE
#error EXODUS_TEST_FILE must be defined
#endif

// -------------------------------------------------------------
// check_preference
// -------------------------------------------------------------

static void
check_preference(Mesh::MeshFactory& mesh_factory, const Mesh::Framework& f)
{
  Mesh::FrameworkPreference pref;
  pref.push_back(f);
  if (Mesh::framework_available(f)) {
    mesh_factory.preference(pref);
    CHECK(mesh_factory.preference().front() == f);
  } else {
    CHECK_THROW(mesh_factory.preference(pref), Mesh::Message);
  }
}


SUITE (MeshFramework)
{

  // This tests setting the Mesh Factory framework preferences.  If
  // only one framework is preferred, and it is not available, an
  // exception should be thrown, while setting preferences
  TEST (PreferenceThrow)
  {
    Epetra_MpiComm comm(MPI_COMM_WORLD);
    bool parallel(comm.NumProc() > 1);
    
    Mesh::MeshFactory mesh_factory(comm);
    Mesh::FrameworkPreference pref(mesh_factory.preference());

    // The Simple framework should always be there
    check_preference(mesh_factory, Mesh::Simple);
    check_preference(mesh_factory, Mesh::MOAB);
    check_preference(mesh_factory, Mesh::STK);
    check_preference(mesh_factory, Mesh::MSTK);
  }
    
  TEST (Generate)
  {
    Epetra_MpiComm comm(MPI_COMM_WORLD);
    bool parallel(comm.NumProc() > 1);
    
    Mesh::FrameworkPreference pref;
    Teuchos::RCP<Mesh_maps_base> mesh;
    Mesh::MeshFactory mesh_factory(comm);

    double x0( 0.0), y0( 0.0), z0( 0.0);
    double x1(10.0), y1(10.0), z1(10.0);
    int nx(10), ny(10), nz(10);

    // The Simple framework is always available, but will only
    // generate in serial

    pref.clear(); pref.push_back(Mesh::Simple);
    mesh_factory.preference(pref);

    if (parallel) {
      CHECK_THROW(mesh = mesh_factory(x0, y0, z0,
                                       x1, y1, z1,
                                       nx, ny, nz),
                  Mesh::Message);
    } else {
      mesh = mesh_factory(x0, y0, z0,
                          x1, y1, z1,
                          nx, ny, nz);
      CHECK(!mesh.is_null());
    }

    // The STK, if available, framework will always generate

    if (framework_available(Mesh::STK)) {
      pref.clear(); pref.push_back(Mesh::STK);
      mesh_factory.preference(pref);
      mesh = mesh_factory(x0, y0, z0,
                          x1, y1, z1,
                          nx, ny, nz);
      CHECK(!mesh.is_null());
    }

    // The MSTK and MOAB frameworks cannot generate

    if (framework_available(Mesh::MSTK)) {
      pref.clear(); pref.push_back(Mesh::MSTK);
      mesh_factory.preference(pref);
      CHECK_THROW(mesh = mesh_factory(x0, y0, z0,
                                      x1, y1, z1,
                                      nx, ny, nz),
                  Mesh::Message);
    }

    if (framework_available(Mesh::MOAB)) {
      pref.clear(); pref.push_back(Mesh::MOAB);
      mesh_factory.preference(pref);
      CHECK_THROW(mesh = mesh_factory(x0, y0, z0,
                                      x1, y1, z1,
                                      nx, ny, nz),
                  Mesh::Message);
    }
  }

  // The Simple framework cannot read anything, even if it exists
  TEST (ReadSimple) {
    Epetra_MpiComm comm(MPI_COMM_WORLD);
    bool parallel(comm.NumProc() > 1);
    
    Mesh::FrameworkPreference pref;
    Teuchos::RCP<Mesh_maps_base> mesh;
    Mesh::MeshFactory mesh_factory(comm);

    pref.clear(); pref.push_back(Mesh::Simple);
    mesh_factory.preference(pref);
    CHECK_THROW(mesh = mesh_factory(BOGUS_TEST_FILE),
                Mesh::Message);
    CHECK_THROW(mesh = mesh_factory(MOAB_TEST_FILE),
                Mesh::Message);
    CHECK_THROW(mesh = mesh_factory(EXODUS_TEST_FILE),
                Mesh::Message);
    CHECK_THROW(mesh = mesh_factory(NEMESIS_TEST_FILE),
                Mesh::Message);
  }

  // Try to read a MOAB HDF5 file, which can only be read by the MOAB
  // framework
  TEST (ReadMOABHDF5)
  {
    Epetra_MpiComm comm(MPI_COMM_WORLD);
    bool parallel(comm.NumProc() > 1);

    Teuchos::RCP<Mesh_maps_base> mesh;
    Mesh::MeshFactory mesh_factory(comm);

    if (Mesh::framework_available(Mesh::MOAB)) {
      mesh = mesh_factory(MOAB_TEST_FILE);
      CHECK(!mesh.is_null());
    } else {
      CHECK_THROW(mesh = mesh_factory(MOAB_TEST_FILE),
                  Mesh::Message);
    }

    // Try it with another framework just for grins

    if (Mesh::framework_available(Mesh::STK)) {
      Mesh::FrameworkPreference pref;
      pref.push_back(Mesh::STK);
      mesh_factory.preference(pref);
      CHECK_THROW(mesh = mesh_factory(MOAB_TEST_FILE),
                  Mesh::Message);
    }
  }

  TEST (ReadExodus) 
  {
    Epetra_MpiComm comm(MPI_COMM_WORLD);
    bool parallel(comm.NumProc() > 1);

    Teuchos::RCP<Mesh_maps_base> mesh;
    Mesh::MeshFactory mesh_factory(comm);

    bool available =
      (Mesh::framework_available(Mesh::STK) && !parallel) ||
      Mesh::framework_available(Mesh::MSTK);

    if (available) {
      mesh = mesh_factory(EXODUS_TEST_FILE);
      CHECK(!mesh.is_null());
    } else {
      CHECK_THROW(mesh = mesh_factory(EXODUS_TEST_FILE),
                  Mesh::Message);
    }
  }
  TEST (ReadNemesis) 
  {
    Epetra_MpiComm comm(MPI_COMM_WORLD);
    bool parallel(comm.NumProc() > 1);
    
    Teuchos::RCP<Mesh_maps_base> mesh;
    Mesh::MeshFactory mesh_factory(comm);
    if (Mesh::framework_available(Mesh::STK) && parallel) {
      mesh = mesh_factory(NEMESIS_TEST_FILE);
      CHECK(!mesh.is_null());
    } else {
      CHECK_THROW(mesh = mesh_factory(NEMESIS_TEST_FILE),
                  Mesh::Message);
    }
  }      

    
}
