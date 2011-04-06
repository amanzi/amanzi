// -------------------------------------------------------------
/**
 * @file   test_mesh_framework.cc
 * @author William A. Perkins
 * @date Tue Mar 22 11:56:46 2011
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created March 14, 2011 by William A. Perkins
// Last Change: Tue Mar 22 11:56:46 2011 by William A. Perkins <d3g096@PE10900.pnl.gov>
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
    
    Mesh::FrameworkPreference pref(Mesh::default_preference());

    CHECK(std::find(pref.begin(), pref.end(), Mesh::Simple) != pref.end());
    
#ifdef HAVE_MOAB_MESH
    CHECK(std::find(pref.begin(), pref.end(), Mesh::MOAB) != pref.end());
#else
    CHECK(std::find(pref.begin(), pref.end(), Mesh::MOAB) == pref.end());
#endif

#ifdef HAVE_STK_MESH
    CHECK(std::find(pref.begin(), pref.end(), Mesh::STK) != pref.end());
#else
    CHECK(std::find(pref.begin(), pref.end(), Mesh::STK) == pref.end());
#endif

#ifdef HAVE_MSTK_MESH
    CHECK(std::find(pref.begin(), pref.end(), Mesh::MSTK) != pref.end());
#else
    CHECK(std::find(pref.begin(), pref.end(), Mesh::MSTK) == pref.end());
#endif

  }

  TEST (AvailablePreference) 
  {
    Mesh::FrameworkPreference pref;
    
    pref.clear();
    pref.push_back(Mesh::MOAB);
    pref = Mesh::available_preference(pref);
#ifdef HAVE_MOAB_MESH
    CHECK(!pref.empty());
#else
    CHECK(pref.empty());
#endif
    
    pref.clear();
    pref.push_back(Mesh::STK);
    pref = Mesh::available_preference(pref);
#ifdef HAVE_STK_MESH
    CHECK(!pref.empty());
#else
    CHECK(pref.empty());
#endif
    
    pref.clear();
    pref.push_back(Mesh::MSTK);
    pref = Mesh::available_preference(pref);
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
    
    CHECK(!Mesh::framework_reads(Mesh::Simple, Mesh::ExodusII, parallel));
    CHECK(!Mesh::framework_reads(Mesh::Simple, Mesh::Nemesis, parallel));
    CHECK(!Mesh::framework_reads(Mesh::Simple, Mesh::MOABHDF5, parallel));

    CHECK(!Mesh::framework_reads(Mesh::MOAB, Mesh::ExodusII, parallel));
    CHECK(!Mesh::framework_reads(Mesh::MOAB, Mesh::Nemesis, parallel));
    CHECK(Mesh::framework_reads(Mesh::MOAB, Mesh::MOABHDF5, parallel));

    if (parallel) {
      CHECK(!Mesh::framework_reads(Mesh::STK, Mesh::ExodusII, parallel));
      CHECK(Mesh::framework_reads(Mesh::STK, Mesh::Nemesis, parallel));
    } else {
      CHECK(Mesh::framework_reads(Mesh::STK, Mesh::ExodusII, parallel));
      CHECK(!Mesh::framework_reads(Mesh::STK, Mesh::Nemesis, parallel));
    }
    CHECK(!Mesh::framework_reads(Mesh::STK, Mesh::MOABHDF5, parallel));

  }

  TEST (Generatability)
  {
    Epetra_MpiComm comm(MPI_COMM_WORLD);
    bool parallel(comm.NumProc() > 1);
    
    CHECK(!Mesh::framework_generates(Mesh::MOAB, parallel));
    CHECK(!Mesh::framework_generates(Mesh::MSTK, parallel));
    CHECK(Mesh::framework_generates(Mesh::STK, parallel));
    if (parallel) {
      CHECK(!Mesh::framework_generates(Mesh::Simple, parallel));
    } 
  }

}
