// -------------------------------------------------------------
/**
 * @file   test_mesh_file.cc
 * @author William A. Perkins
 * @date Mon Mar 14 10:24:48 2011
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
// Last Change: Mon Mar 14 10:24:48 2011 by William A. Perkins <d3g096@PE10900.pnl.gov>
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

    Mesh::Format f;
    try {
      f = Mesh::file_format(comm, fname);
    } catch (const Mesh::Message& e) {
      throw e;
    }

    CHECK(f == Mesh::ExodusII);
  }

  TEST (Nemesis) 
  {
    Epetra_MpiComm comm(MPI_COMM_WORLD);

    // NEMESIS_TEST_FILE is macro defined by cmake
    std::string fname(NEMESIS_TEST_FILE); 
    
    Mesh::Format f;
    if (comm.NumProc() > 1 && comm.NumProc() <= 4) {
      int ierr[1];
      ierr[0] = 0;
      try {
        f = Mesh::file_format(comm, fname);
      } catch (const Mesh::Message& e) {
        throw e;
      }

      CHECK(f == Mesh::Nemesis);
    } else {
      CHECK_THROW(f = Mesh::file_format(comm, fname), Mesh::FileMessage);
    }  
  }
   
  TEST (MOABHD5) 
  {
    Epetra_MpiComm comm(MPI_COMM_WORLD);

    // MOAB_TEST_FILE is macro defined by cmake
    std::string fname(MOAB_TEST_FILE); 

    Mesh::Format f;
    try {
      f = Mesh::file_format(comm, fname);
    } catch (const Mesh::Message& e) {
      throw e;
    }

    CHECK(f == Mesh::MOABHDF5);
  }

  TEST (PathFailure) 
  {
    Epetra_MpiComm comm(MPI_COMM_WORLD);

    std::string fname("/some/bogus/path.exo"); 

    Mesh::Format f;

    CHECK_THROW(Mesh::file_format(comm, fname), Mesh::FileMessage);
  }    

  TEST (MagicNumberFailure)
  {
    Epetra_MpiComm comm(MPI_COMM_WORLD);

    std::string fname(BOGUS_TEST_FILE); 

    Mesh::Format f;

    CHECK_THROW(Mesh::file_format(comm, fname), Mesh::FileMessage);
  }
    
}
    
