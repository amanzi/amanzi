// -------------------------------------------------------------
/**
 * @file   verify_stk_mesh.cc
 * @author William A. Perkins
 * @date Mon Dec 13 10:35:09 2010
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created December 13, 2010 by William A. Perkins
// Last Change: Mon Dec 13 10:35:09 2010 by William A. Perkins <d3g096@PE10900.pnl.gov>
// -------------------------------------------------------------


#include <mpi.h>
#include <iostream>

#include "Teuchos_GlobalMPISession.hpp"
#include "Epetra_MpiComm.h"
#include "Epetra_SerialComm.h"

#include "Exodus_readers.hh"
#include "Parallel_Exodus_file.hh"
#include "Mesh_factory.hh"
#include "Mesh_maps_stk.hh"

#include "MeshAudit.hh"

int main (int argc, char* argv[])
{

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  Epetra_MpiComm comm(MPI_COMM_WORLD);
  const int nproc(comm.NumProc());
  const int me(comm.MyPID());

  if (argc < 2) {
    std::cerr << "Usage: " << argv[0] << " name" << std::endl;
    return 3;
  }

  // the first, and only, command line argument is a file name; when
  // run in parallel this is a base name to which the number of
  // processes and process id is added, e.g., basename.N.P; if run
  // serially the file name is the complete name of the file

  std::string filename(argv[1]);

  Teuchos::RCP<Mesh_maps_base> 
    maps(new STK_mesh::Mesh_maps_stk(comm, filename.c_str()));

  int status;

  if (nproc == 1) {
    MeshAudit audit(maps);
    status = audit.Verify();
  } else {
    std::ostringstream ofile;
    ofile << "mesh_audit_" << std::setfill('0') << std::setw(4) << me << ".txt";
    std::ofstream ofs(ofile.str().c_str());
    if (me == 0)
      std::cout << "Writing results to " << ofile.str() << ", etc." << std::endl;
    MeshAudit audit(maps, ofs);
    status = audit.Verify();
  }

  if (me == 0) {
    std::cout << filename << ": " << (status ? "has errors" : "OK") << std::endl;
  }

  return status;
}
