// -------------------------------------------------------------
/**
 * @file   verify_stk_mesh.cc
 * @author William A. Perkins
 * @date Mon Dec 13 08:02:50 2010
 * 
 * @brief  
 * 
 * 
 */
// -------------------------------------------------------------
// -------------------------------------------------------------
// Created December 13, 2010 by William A. Perkins
// Last Change: Mon Dec 13 08:02:50 2010 by William A. Perkins <d3g096@PE10900.pnl.gov>
// -------------------------------------------------------------


#include <mpi.h>
#include <iostream>

#include "Teuchos_GlobalMPISession.hpp"
#include "Epetra_MpiComm.h"
#include "Epetra_SerialComm.h"

#include "Exodus_Readers.hh"
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

  Teuchos::RCP<Mesh_data::Data> meshdata;
  STK_mesh::Mesh_p mesh;
  STK_mesh::Mesh_factory mf(comm.Comm(), 1000);
  Mesh_data::Fields nofields;
  Teuchos::RCP<Mesh_maps_base> maps;

  if (nproc == 1) {
    meshdata.reset(ExodusII::read_exodus_file(filename.c_str()));
    mesh.reset(mf.build_mesh(*meshdata, nofields));
    mesh->summary(std::cout);
    maps.reset(new STK_mesh::Mesh_maps_stk(mesh));

    MeshAudit audit(maps);
    return audit.Verify();
  } else {
    ExodusII::Parallel_Exodus_file thefile(comm, filename.c_str());
    meshdata = thefile.read_mesh();
    Teuchos::RCP<Epetra_Map> cmap(thefile.cellmap());
    Teuchos::RCP<Epetra_Map> vmap(thefile.vertexmap());

    mesh.reset(mf.build_mesh(*meshdata, *cmap, *vmap, nofields));
    mesh->summary(std::cout);
    maps.reset(new STK_mesh::Mesh_maps_stk(mesh));

    std::ostringstream ofile;
    ofile << "mesh_audit_" << std::setfill('0') << std::setw(4) << me << ".txt";
    std::ofstream ofs(ofile.str().c_str());
    if (me == 0)
      std::cout << "Writing results to " << ofile.str() << ", etc." << std::endl;
    MeshAudit audit(maps, ofs);
    return audit.Verify();
  }

  return 0;
}
