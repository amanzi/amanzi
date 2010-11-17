// -------------------------------------------------------------
/**
 * @file   test_pread_2.cc
 * @author William A. Perkins
 * @date Wed Nov 17 09:04:28 2010
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created November 15, 2010 by William A. Perkins
// Last Change: Wed Nov 17 09:04:28 2010 by William A. Perkins <d3g096@PE10900.pnl.gov>
// -------------------------------------------------------------

#include <iostream>
#include <UnitTest++.h>

#include <Epetra_Map.h>
#include <Epetra_MpiComm.h>

#include "dbc.hh"
#include "../Parallel_Exodus_file.hh"

extern std::string split_file_path(const std::string& fname);
extern void checkit(ExodusII::Parallel_Exodus_file & thefile);

SUITE (Exodus_2_Proc)
{
  TEST (hex_4x4x4_ss)
  {
    std::string bname(split_file_path("hex_4x4x4_ss.par").c_str());
    
    Epetra_MpiComm comm(MPI_COMM_WORLD);

    CHECK_EQUAL(comm.NumProc(), 2);
    
    ExodusII::Parallel_Exodus_file thefile(comm, bname);
    checkit(thefile);
  }


  TEST (htc_rad_test_random)
  {
    std::string bname(split_file_path("htc_rad_test-random.par").c_str());
    
    Epetra_MpiComm comm(MPI_COMM_WORLD);

    CHECK_EQUAL(comm.NumProc(), 2);
    
    ExodusII::Parallel_Exodus_file thefile(comm, bname);
    checkit(thefile);
      
  }

  TEST (hex_11x11x11_ss)
  {
    std::string bname(split_file_path("hex_11x11x11_ss.par").c_str());
    
    Epetra_MpiComm comm(MPI_COMM_WORLD);

    CHECK_EQUAL(comm.NumProc(), 2);
    
    ExodusII::Parallel_Exodus_file thefile(comm, bname);
    checkit(thefile);
      
  }

}
