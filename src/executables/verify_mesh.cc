/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
// -------------------------------------------------------------
/**
 * @file   verify_mesh.cc
 * @author William A. Perkins
 * @date Mon Jun 20 11:37:31 2011
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created December 13, 2010 by William A. Perkins
// Last Change: Mon Jun 20 11:37:31 2011 by William A. Perkins <d3g096@PE10900.pnl.gov>
// -------------------------------------------------------------


#include <mpi.h>
#include <iostream>

#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
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

  // handle command line

  Teuchos::CommandLineProcessor CLP;
  
  CLP.setDocString("reads mesh file or file set and does a series of checks\n");

  const Amanzi::AmanziMesh::Framework frameworks[] = {  
    Amanzi::AmanziMesh::MOAB, 
    Amanzi::AmanziMesh::STKMESH, 
    Amanzi::AmanziMesh::MSTK 
  };
  const char *framework_names[] = {
    "MOAB", "stk::mesh", "MSTK"
  };

  const int numframeworks = sizeof(frameworks)/sizeof(Amanzi::AmanziMesh::Framework);
  
  Amanzi::AmanziMesh::Framework the_framework(Amanzi::AmanziMesh::MOAB);
  CLP.setOption("framework", &the_framework,
                numframeworks, frameworks, framework_names,
                "mesh framework preference", false);

  std::string filename;
  CLP.setOption("file", &filename, "name of file or file set", true);

  CLP.throwExceptions(false);

  int ierr(0);
  Teuchos::CommandLineProcessor::EParseCommandLineReturn parseReturn;
  try {
    parseReturn = CLP.parse(argc, argv);
  } catch (const std::exception &e) {
    std::cerr << "error: " << e.what() << std::endl;
    ierr++;
  }

  comm.SumAll(&ierr, &ierr, 1);

  if (ierr > 0) {
    return 1;
  }
  
  if (parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) {
    return 0;
  }
  
  if (parseReturn != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
    return 1;
  }

  // the first, and only, command line argument is a file name. Three
  // types are supported depending on which frameworks are compiled in

  Amanzi::AmanziMesh::MeshFactory factory(comm);
  Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh;
  
  ierr = 0;
  try {
    Amanzi::AmanziMesh::FrameworkPreference prefs(factory.preference());
    if (the_framework !=  Amanzi::AmanziMesh::Simple) {
      prefs.clear(); 
      prefs.push_back(the_framework);
    } 

    if (me == 0) {
      std::cerr << "Attempting to read \"" << filename << "\" with ";
      for (Amanzi::AmanziMesh::FrameworkPreference::iterator i = prefs.begin();
           i != prefs.end(); i++) {
        std::cerr << Amanzi::AmanziMesh::framework_name(*i) << ", ";
      }
      std::cerr << std::endl;
    }
    factory.preference(prefs);
    mesh = factory(filename);
  } catch (const Amanzi::AmanziMesh::Message& e) {
    std::cerr << argv[0] << ": mesh error: " << e.what() << std::endl;
    ierr++;
  } catch (const std::exception& e) {
    std::cerr << argv[0] << ": error: " << e.what() << std::endl;
    ierr++;
  }

  comm.SumAll(&ierr, &ierr, 1);

  if (ierr > 0) {
    return 3;
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
