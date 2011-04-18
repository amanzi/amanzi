#include "Mesh_maps_moab.hh"

#include <mpi.h>
#include <iostream>

#include "Teuchos_GlobalMPISession.hpp"
#include "Epetra_MpiComm.h"
#include "Epetra_SerialComm.h"
#include "Mesh_maps_base.hh"

#include "MeshAudit.hh"

using namespace std;

int main (int argc, char* argv[])
{
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  Epetra_MpiComm comm(MPI_COMM_WORLD);
  
  // The program takes a single argument that is the name of the mesh file.
  
  if (argc != 2) {
    cerr << "Usage: " << argv[0] << "file" << endl;
    return 2;
  }
  
  string filename(argv[1]);

  Teuchos::RCP<Mesh_maps_base> mesh(new Mesh_maps_moab(filename.c_str(), MPI_COMM_WORLD));
  
  int status;

  if (comm.NumProc() == 1) {
    MeshAudit audit(mesh);
    status = audit.Verify();
  } else {
    ostringstream ofile;
    ofile << "mesh_audit_" << setfill('0') << setw(4) << comm.MyPID() << ".txt";
    ofstream ofs(ofile.str().c_str());
    if (comm.MyPID() == 0)
      cout << "Writing results to " << ofile.str() << ", etc." << endl;
    MeshAudit audit(mesh, ofs);
    status = audit.Verify();
  }
  
  if (comm.MyPID() == 0) {
    cerr << filename << ": " << (status ? "has errors" : "OK") << endl;
  }
  
  return status;
}
