// -------------------------------------------------------------
/**
 * @file   test_mesh_file.cc
 * @author William A. Perkins
 * @date Mon May 16 14:03:23 2011
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
// Battelle Memorial Institute
// Pacific Northwest Laboratory
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created March 14, 2011 by William A. Perkins
// Last Change: Mon May 16 14:03:23 2011 by William A. Perkins <d3g096@PE10900.pnl.gov>
// -------------------------------------------------------------


static const char* SCCS_ID = "$Id$ Battelle PNL";

#include <iostream>
#include <UnitTest++.h>

#include <Epetra_MpiComm.h>

#include "dbc.hh"
#include "../MeshFileType.hh"
#include "../MeshException.hh"

SUITE (MeshFileType)
{
  TEST (ExodusII)
  {
    Epetra_MpiComm comm(MPI_COMM_WORLD);

    // EXODUS_TEST_FILE is macro defined by cmake
    std::string fname(EXODUS_TEST_FILE); 

    Amanzi::AmanziMesh::Format f;
    try {
      f = Amanzi::AmanziMesh::file_format(comm, fname);
    } catch (const Amanzi::AmanziMesh::Message& e) {
      throw e;
    }

    CHECK(f == Amanzi::AmanziMesh::ExodusII);
  }

  TEST (Nemesis) 
  {
    Epetra_MpiComm comm(MPI_COMM_WORLD);

    // NEMESIS_TEST_FILE is macro defined by cmake
    std::string fname(NEMESIS_TEST_FILE); 
    
    Amanzi::AmanziMesh::Format f;
    if (comm.NumProc() > 1 && comm.NumProc() <= 4) {
      int ierr[1];
      ierr[0] = 0;
      try {
        f = Amanzi::AmanziMesh::file_format(comm, fname);
      } catch (const Amanzi::AmanziMesh::Message& e) {
        throw e;
      }

      CHECK(f == Amanzi::AmanziMesh::Nemesis);
    } else {
      CHECK_THROW(f = Amanzi::AmanziMesh::file_format(comm, fname), 
                  Amanzi::AmanziMesh::FileMessage);
    }  
  }
   
  TEST (MOABHD5) 
  {
    Epetra_MpiComm comm(MPI_COMM_WORLD);

    // MOAB_TEST_FILE is macro defined by cmake
    std::string fname(MOAB_TEST_FILE); 

    Amanzi::AmanziMesh::Format f;
    try {
      f = Amanzi::AmanziMesh::file_format(comm, fname);
    } catch (const Amanzi::AmanziMesh::Message& e) {
      throw e;
    }

    CHECK(f == Amanzi::AmanziMesh::MOABHDF5);
  }

  TEST (PathFailure) 
  {
    Epetra_MpiComm comm(MPI_COMM_WORLD);

    std::string fname("/some/bogus/path.exo"); 

    Amanzi::AmanziMesh::Format f;

    CHECK_THROW(Amanzi::AmanziMesh::file_format(comm, fname), 
                Amanzi::AmanziMesh::FileMessage);
  }    

  TEST (MagicNumberFailure)
  {
    Epetra_MpiComm comm(MPI_COMM_WORLD);

    std::string fname(BOGUS_TEST_FILE); 

    Amanzi::AmanziMesh::Format f;

    CHECK_THROW(Amanzi::AmanziMesh::file_format(comm, fname), 
                Amanzi::AmanziMesh::FileMessage);
  }
    
}
    
