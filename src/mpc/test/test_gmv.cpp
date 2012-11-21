#include <iostream>
#include "stdlib.h"
#include "math.h"
#include "UnitTest++.h"
#include "MeshFactory.hh"
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
  
  Amanzi::AmanziMesh::FrameworkPreference pref;
  pref.clear();
  pref.push_back(Amanzi::AmanziMesh::Simple);

  Amanzi::AmanziMesh::MeshFactory meshfactory(comm);
  meshfactory.preference(pref);
  Teuchos::RCP<Amanzi::AmanziMesh::Mesh> MMS =  meshfactory(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 3, 4, 2);

  State S(1,MMS);

  std::string gmv_filename = "test_gmv.gmv";
  
  S.write_gmv(gmv_filename);
      
}


