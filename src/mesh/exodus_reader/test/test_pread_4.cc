// -------------------------------------------------------------
/**
 * @file   test_pread_2.cc
 * @author William A. Perkins
 * @date Mon May  2 10:18:48 2011
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created November 15, 2010 by William A. Perkins
// Last Change: Mon May  2 10:18:48 2011 by William A. Perkins <d3g096@PE10900.pnl.gov>
// -------------------------------------------------------------

#include <iostream>
#include <UnitTest++.h>

#include <Epetra_MpiComm.h>

#include "dbc.hh"
#include "../Parallel_Exodus_file.hh"

extern std::string split_file_path(const std::string& fname);
extern void checkit(Amanzi::Exodus::Parallel_Exodus_file & thefile);

SUITE (Exodus_4_Proc)
{
  TEST (htc_rad_test_random)
  {
    std::string bname(split_file_path("htc_rad_test-random.par").c_str());
    
    Epetra_MpiComm comm(MPI_COMM_WORLD);

    CHECK_EQUAL(comm.NumProc(), 4);
    
    Amanzi::Exodus::Parallel_Exodus_file thefile(comm, bname);
    checkit(thefile);
  }

  TEST (hex_10x10x10_ss)
  {
    std::string bname(split_file_path("hex_10x10x10_ss.par").c_str());
    
    Epetra_MpiComm comm(MPI_COMM_WORLD);

    CHECK_EQUAL(comm.NumProc(), 4);
    
    Amanzi::Exodus::Parallel_Exodus_file thefile(comm, bname);
    checkit(thefile);
  }

  TEST (twoblktet_ss)
  {
    std::string bname(split_file_path("twoblktet_ss.par").c_str());
    
    Epetra_MpiComm comm(MPI_COMM_WORLD);

    CHECK_EQUAL(comm.NumProc(), 4);
    
    Amanzi::Exodus::Parallel_Exodus_file thefile(comm, bname);
    checkit(thefile);
  }

  TEST (prism)
  {
    std::string bname(split_file_path("prism.par").c_str());
    
    Epetra_MpiComm comm(MPI_COMM_WORLD);

    CHECK_EQUAL(comm.NumProc(), 4);
    
    Amanzi::Exodus::Parallel_Exodus_file thefile(comm, bname);
    checkit(thefile);
      
  }

  TEST (mixed_coarse)
  {
    std::string bname(split_file_path("mixed-coarse.par").c_str());
    
    Epetra_MpiComm comm(MPI_COMM_WORLD);

    CHECK_EQUAL(comm.NumProc(), 4);
    
    Amanzi::Exodus::Parallel_Exodus_file thefile(comm, bname);
    checkit(thefile);
      
  }

  TEST (mixed)
  {
    std::string bname(split_file_path("mixed.par").c_str());
    
    Epetra_MpiComm comm(MPI_COMM_WORLD);

    CHECK_EQUAL(comm.NumProc(), 4);
    
    Amanzi::Exodus::Parallel_Exodus_file thefile(comm, bname);
    checkit(thefile);
      
  }
}
