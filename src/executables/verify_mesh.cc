/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
// -------------------------------------------------------------
/**
 * @file   verify_mesh.cc
 * @author William A. Perkins
 * @date Wed May 18 12:33:22 2011
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created December 13, 2010 by William A. Perkins
// Last Change: Wed May 18 12:33:22 2011 by William A. Perkins <d3g096@PE10900.pnl.gov>
// -------------------------------------------------------------


#include <mpi.h>
#include <iostream>

#include "Teuchos_GlobalMPISession.hpp"
#include "Epetra_MpiComm.h"
#include "Epetra_SerialComm.h"

#include "MeshFactory.hh"
#include "MeshAudit.hh"
#include "MeshException.hh"

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

  // the first, and only, command line argument is a file name. Three
  // types are supported depending on which frameworks are compiled in

  std::string filename(argv[1]);

  Amanzi::AmanziMesh::MeshFactory factory(comm);
  Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh;

  try {
    mesh = factory(filename);
  } catch (const Amanzi::AmanziMesh::Message& e) {
    std::cerr << argv[0] << ": mesh error: " << e.what() << std::endl;
    exit(3);
  } catch (const std::exception& e) {
    std::cerr << argv[0] << ": error: " << e.what() << std::endl;
    exit(3);
  }

  int status;

  if (nproc == 1) {
    Amanzi::MeshAudit audit(mesh);
    status = audit.Verify();
  } else {
    std::ostringstream ofile;
    ofile << "mesh_audit_" << std::setfill('0') << std::setw(4) << me << ".txt";
    std::ofstream ofs(ofile.str().c_str());
    if (me == 0)
      std::cout << "Writing results to " << ofile.str() << ", etc." << std::endl;
    Amanzi::MeshAudit audit(mesh, ofs);
    status = audit.Verify();
  }

  if (me == 0) {
    std::cout << filename << ": " << (status ? "has errors" : "OK") << std::endl;
  }

  return status;
}
