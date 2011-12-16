// -------------------------------------------------------------
/**
 * @file   test_mesh_framework.cc
 * @author William A. Perkins
 * @date Tue Oct  4 06:15:36 2011
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created March 14, 2011 by William A. Perkins
// Last Change: Tue Oct  4 06:15:36 2011 by William A. Perkins <d3g096@PE10900.pnl.gov>
// -------------------------------------------------------------

#include <iostream>
#include <UnitTest++.h>

#include <Epetra_MpiComm.h>

#include "dbc.hh"
#include "../MeshFramework.hh"
#include "../FrameworkTraits.hh"

SUITE (Framework) 
{
  TEST (DefaultPreference) {
    
    Amanzi::AmanziMesh::FrameworkPreference pref(Amanzi::AmanziMesh::default_preference());

    CHECK(std::find(pref.begin(), pref.end(), Amanzi::AmanziMesh::Simple) != pref.end());
    
#ifdef HAVE_MOAB_MESH
    CHECK(std::find(pref.begin(), pref.end(), Amanzi::AmanziMesh::MOAB) != pref.end());
#else
    CHECK(std::find(pref.begin(), pref.end(), Amanzi::AmanziMesh::MOAB) == pref.end());
#endif

#ifdef HAVE_STK_MESH
    CHECK(std::find(pref.begin(), pref.end(), Amanzi::AmanziMesh::STKMESH) != pref.end());
#else
    CHECK(std::find(pref.begin(), pref.end(), Amanzi::AmanziMesh::STKMESH) == pref.end());
#endif

#ifdef HAVE_MSTK_MESH
    CHECK(std::find(pref.begin(), pref.end(), Amanzi::AmanziMesh::MSTK) != pref.end());
#else
    CHECK(std::find(pref.begin(), pref.end(), Amanzi::AmanziMesh::MSTK) == pref.end());
#endif

  }

  TEST (AvailablePreference) 
  {
    Amanzi::AmanziMesh::FrameworkPreference pref;
    
    pref.clear();
    pref.push_back(Amanzi::AmanziMesh::MOAB);
    pref = Amanzi::AmanziMesh::available_preference(pref);
#ifdef HAVE_MOAB_MESH
    CHECK(!pref.empty());
#else
    CHECK(pref.empty());
#endif
    
    pref.clear();
    pref.push_back(Amanzi::AmanziMesh::STKMESH);
    pref = Amanzi::AmanziMesh::available_preference(pref);
#ifdef HAVE_STK_MESH
    CHECK(!pref.empty());
#else
    CHECK(pref.empty());
#endif
    
    pref.clear();
    pref.push_back(Amanzi::AmanziMesh::MSTK);
    pref = Amanzi::AmanziMesh::available_preference(pref);
#ifdef HAVE_MSTK_MESH
    CHECK(!pref.empty());
#else
    CHECK(pref.empty());
#endif
    
  }

  TEST (Readability)
  {
    Epetra_MpiComm comm(MPI_COMM_WORLD);
    bool parallel(comm.NumProc() > 1);
    
    CHECK(!Amanzi::AmanziMesh::framework_reads(Amanzi::AmanziMesh::Simple, Amanzi::AmanziMesh::ExodusII, parallel));
    CHECK(!Amanzi::AmanziMesh::framework_reads(Amanzi::AmanziMesh::Simple, Amanzi::AmanziMesh::Nemesis, parallel));
    CHECK(!Amanzi::AmanziMesh::framework_reads(Amanzi::AmanziMesh::Simple, Amanzi::AmanziMesh::MOABHDF5, parallel));

    CHECK(Amanzi::AmanziMesh::framework_reads(Amanzi::AmanziMesh::MOAB, Amanzi::AmanziMesh::ExodusII, parallel));
    CHECK(!Amanzi::AmanziMesh::framework_reads(Amanzi::AmanziMesh::MOAB, Amanzi::AmanziMesh::Nemesis, parallel));
    CHECK(Amanzi::AmanziMesh::framework_reads(Amanzi::AmanziMesh::MOAB, Amanzi::AmanziMesh::MOABHDF5, parallel));

    if (parallel) {
      CHECK(!Amanzi::AmanziMesh::framework_reads(Amanzi::AmanziMesh::STKMESH, Amanzi::AmanziMesh::ExodusII, parallel));
      CHECK(Amanzi::AmanziMesh::framework_reads(Amanzi::AmanziMesh::STKMESH, Amanzi::AmanziMesh::Nemesis, parallel));
    } else {
      CHECK(Amanzi::AmanziMesh::framework_reads(Amanzi::AmanziMesh::STKMESH, Amanzi::AmanziMesh::ExodusII, parallel));
      CHECK(!Amanzi::AmanziMesh::framework_reads(Amanzi::AmanziMesh::STKMESH, Amanzi::AmanziMesh::Nemesis, parallel));
    }
    CHECK(!Amanzi::AmanziMesh::framework_reads(Amanzi::AmanziMesh::STKMESH, Amanzi::AmanziMesh::MOABHDF5, parallel));

  }

  TEST (Generatability)
  {
    Epetra_MpiComm comm(MPI_COMM_WORLD);
    bool parallel(comm.NumProc() > 1);
    
    CHECK(!Amanzi::AmanziMesh::framework_generates(Amanzi::AmanziMesh::MOAB, parallel));
    CHECK(Amanzi::AmanziMesh::framework_generates(Amanzi::AmanziMesh::MSTK, parallel));
    CHECK(Amanzi::AmanziMesh::framework_generates(Amanzi::AmanziMesh::STKMESH, parallel));
    if (parallel) {
      CHECK(!Amanzi::AmanziMesh::framework_generates(Amanzi::AmanziMesh::Simple, parallel));
    } 
  }

}
