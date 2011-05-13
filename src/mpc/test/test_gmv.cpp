#include <iostream>
#include "stdlib.h"
#include "math.h"
#include "UnitTest++.h"
#include "Mesh_simple.hh"
#include <Epetra_Comm.h>
#include <Epetra_MpiComm.h>
#include "Epetra_SerialComm.h"
#include "../State.hpp"


TEST(GMV) {

#ifdef HAVE_MPI
  Epetra_MpiComm *comm = new Epetra_MpiComm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm *comm = new Epetra_SerialComm();
#endif
  
  Teuchos::RCP<Mesh_simple> MMS = 
    Teuchos::rcp(new Mesh_simple(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 3, 4, 2, comm ));

  State S(1,MMS);

  std::string gmv_filename = "test_gmv.gmv";
  
  S.write_gmv(gmv_filename);
      
}


