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

#include "AmanziTypes.hh"
#include "Teuchos_DefaultMpiComm.hpp"

#include "dbc.hh"
#include "../MeshFramework.hh"
#include "../FrameworkTraits.hh"

using namespace Amanzi;

SUITE (Framework) 
{
  TEST (DefaultPreference) {
    
    AmanziMesh::FrameworkPreference pref(AmanziMesh::default_preference());

    CHECK(std::find(pref.begin(), pref.end(), AmanziMesh::Simple) != pref.end());
    
#ifdef HAVE_MOAB_MESH
    CHECK(std::find(pref.begin(), pref.end(), AmanziMesh::MOAB) != pref.end());
#else
    CHECK(std::find(pref.begin(), pref.end(), AmanziMesh::MOAB) == pref.end());
#endif

#ifdef HAVE_STK_MESH
    CHECK(std::find(pref.begin(), pref.end(), AmanziMesh::STKMESH) != pref.end());
#else
    CHECK(std::find(pref.begin(), pref.end(), AmanziMesh::STKMESH) == pref.end());
#endif

#ifdef HAVE_MSTK_MESH
    CHECK(std::find(pref.begin(), pref.end(), AmanziMesh::MSTK) != pref.end());
#else
    CHECK(std::find(pref.begin(), pref.end(), AmanziMesh::MSTK) == pref.end());
#endif

  }

  TEST (AvailablePreference) 
  {
    AmanziMesh::FrameworkPreference pref;
    
    pref.clear();
    pref.push_back(AmanziMesh::MOAB);
    pref = AmanziMesh::available_preference(pref);
#ifdef HAVE_MOAB_MESH
    CHECK(!pref.empty());
#else
    CHECK(pref.empty());
#endif
    
    pref.clear();
    pref.push_back(AmanziMesh::STKMESH);
    pref = AmanziMesh::available_preference(pref);
#ifdef HAVE_STK_MESH
    CHECK(!pref.empty());
#else
    CHECK(pref.empty());
#endif
    
    pref.clear();
    pref.push_back(AmanziMesh::MSTK);
    pref = AmanziMesh::available_preference(pref);
#ifdef HAVE_MSTK_MESH
    CHECK(!pref.empty());
#else
    CHECK(pref.empty());
#endif
    
  }

  TEST (Readability)
  {
    auto comm = Comm_ptr_type(new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
    bool parallel(comm->getSize() > 1);
    
    CHECK(!AmanziMesh::framework_reads(AmanziMesh::Simple, AmanziMesh::ExodusII, parallel));
    CHECK(!AmanziMesh::framework_reads(AmanziMesh::Simple, AmanziMesh::Nemesis, parallel));
    CHECK(!AmanziMesh::framework_reads(AmanziMesh::Simple, AmanziMesh::MOABHDF5, parallel));

    CHECK(AmanziMesh::framework_reads(AmanziMesh::MOAB, AmanziMesh::ExodusII, parallel));
    CHECK(!AmanziMesh::framework_reads(AmanziMesh::MOAB, AmanziMesh::Nemesis, parallel));
    CHECK(AmanziMesh::framework_reads(AmanziMesh::MOAB, AmanziMesh::MOABHDF5, parallel));

    if (parallel) {
      CHECK(!AmanziMesh::framework_reads(AmanziMesh::STKMESH, AmanziMesh::ExodusII, parallel));
      CHECK(AmanziMesh::framework_reads(AmanziMesh::STKMESH, AmanziMesh::Nemesis, parallel));
    } else {
      CHECK(AmanziMesh::framework_reads(AmanziMesh::STKMESH, AmanziMesh::ExodusII, parallel));
      CHECK(!AmanziMesh::framework_reads(AmanziMesh::STKMESH, AmanziMesh::Nemesis, parallel));
    }
    CHECK(!AmanziMesh::framework_reads(AmanziMesh::STKMESH, AmanziMesh::MOABHDF5, parallel));

  }

  TEST (Generatability)
  {
    auto comm = Comm_ptr_type(new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
    bool parallel(comm->getSize() > 1);
    
    CHECK(!AmanziMesh::framework_generates(AmanziMesh::MOAB, parallel,3));
    CHECK(AmanziMesh::framework_generates(AmanziMesh::MSTK, parallel,3));
    CHECK(AmanziMesh::framework_generates(AmanziMesh::STKMESH, parallel,3));
    if (parallel) {
      CHECK(!AmanziMesh::framework_generates(AmanziMesh::Simple, parallel,3));
    } 


    //    CHECK(!AmanziMesh::framework_generates(AmanziMesh::MOAB, parallel,2));
    //    CHECK(AmanziMesh::framework_generates(AmanziMesh::MSTK, parallel,2));
    //    CHECK(!AmanziMesh::framework_generates(AmanziMesh::STKMESH, parallel,2));
    //    if (parallel) {
    //      CHECK(!AmanziMesh::framework_generates(AmanziMesh::Simple, parallel,2));
    //    } 

  }

}
