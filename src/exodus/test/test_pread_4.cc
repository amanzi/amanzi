// -------------------------------------------------------------
/**
 * @file   test_pread_2.cc
 * @author William A. Perkins
 * @date Tue Nov 16 11:28:40 2010
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created November 15, 2010 by William A. Perkins
// Last Change: Tue Nov 16 11:28:40 2010 by William A. Perkins <d3g096@PE10900.pnl.gov>
// -------------------------------------------------------------

#include <iostream>
#include <UnitTest++.h>

#include <Epetra_MpiComm.h>

#include "dbc.hh"
#include "../Parallel_Exodus_file.hh"

extern void checkit(ExodusII::Parallel_Exodus_file & thefile);

SUITE (Exodus_4_Proc)
{
  TEST (htc_rad_test_random)
  {
    std::string bname("test_files/split1/htc_rad_test-random.par");
    
    Epetra_MpiComm comm(MPI_COMM_WORLD);

    CHECK_EQUAL(comm.NumProc(), 4);
    
    ExodusII::Parallel_Exodus_file thefile(comm, bname);
    checkit(thefile);
  }

  TEST (hex_11x11x11_ss)
  {
    std::string bname("test_files/split1/hex_11x11x11_ss.par");
    
    Epetra_MpiComm comm(MPI_COMM_WORLD);

    CHECK_EQUAL(comm.NumProc(), 4);
    
    ExodusII::Parallel_Exodus_file thefile(comm, bname);
    checkit(thefile);
  }

  TEST (twoblktet_ss)
  {
    std::string bname("test_files/split1/twoblktet_ss.par");
    
    Epetra_MpiComm comm(MPI_COMM_WORLD);

    CHECK_EQUAL(comm.NumProc(), 4);
    
    ExodusII::Parallel_Exodus_file thefile(comm, bname);
    checkit(thefile);
  }

}
