#include <iostream>
#include "stdlib.h"
#include "math.h"
#include "UnitTest++.h"
#include "MeshFactory.hh"
#include <Epetra_Comm.h>
#include <Epetra_MpiComm.h>
#include "Epetra_SerialComm.h"
#include "State.hh"


TEST(GMV) {

#ifdef HAVE_MPI
  auto comm = Amanzi::getDefaultComm();
#else
  Epetra_SerialComm *comm = new Epetra_SerialComm();
#endif
  
  Amanzi::AmanziMesh::Preference pref;
  pref.clear();
  pref.push_back(Amanzi::AmanziMesh::Framework::SIMPLE);

  Amanzi::AmanziMesh::MeshFactory meshfactory(comm);
  meshfactory.set_preference(pref);
  Teuchos::RCP<Amanzi::AmanziMesh::Mesh> MMS =  meshfactory.create(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 3, 4, 2);

  State S(1,MMS);

  std::string gmv_filename = "test_gmv.gmv";
  
  S.write_gmv(gmv_filename);
      
}


