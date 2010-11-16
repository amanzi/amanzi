#include "Mesh_maps_moab.hh"

#include <mpi.h>
#include <iostream>

#include "Teuchos_GlobalMPISession.hpp"
#include "Epetra_MpiComm.h"
#include "Epetra_SerialComm.h"
#include "Mesh_maps_base.hh"

#include "MeshAudit.hh"

//using namespace std;

int main (int argc, char* argv[])
{

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

#ifdef HAVE_MPI
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif

  Teuchos::RCP<Mesh_maps_base> mesh(new Mesh_maps_moab(argv[1], MPI_COMM_WORLD));

  if (comm.NumProc() == 1) {
    MeshAudit audit(mesh);
    return audit.Verify();
  } else {
    std::ostringstream ofile;
    ofile << "mesh_audit_" << std::setfill('0') << std::setw(4) << comm.MyPID() << ".txt";
    std::ofstream ofs(ofile.str().c_str());
    if (comm.MyPID() == 0)
      std::cout << "Writing results to " << ofile.str() << ", etc." << std::endl;
    MeshAudit audit(mesh, ofs);
    return audit.Verify();
  }
}
