#include <cstdlib>
#include <iostream>

#include "UnitTest++.h"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "Mesh_MSTK.hh"

#include "../deform.hh"

/* **************************************************************** */
TEST(ADVANCE_WITH_SIMPLE) {
  using namespace Amanzi;
  using namespace Amanzi::AmanziMesh;
  using namespace Amanzi::AmanziGeometry;

  // first message
  std::cout << "Test: DeformMesh using a mesh from Exodus file" << endl;

  // sequential/parallel communicator object
#ifdef HAVE_MPI
  Epetra_MpiComm *comm = new Epetra_MpiComm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm *comm = new Epetra_SerialComm();
#endif
  
  // read parameter list from input XML file
  Teuchos::ParameterList plist;
  string xml_fname = "test/simple_mesh_deform_test_2D.xml";
  Teuchos::Ptr<Teuchos::ParameterList> plist_p(&plist);
  updateParametersFromXmlFile(xml_fname, plist_p);

  // intantiate the mesh: get the region list
  Teuchos::ParameterList region_list = 
    plist.get<Teuchos::ParameterList>("Regions");
  
  // intantiate the mesh: build the geometric model
  const int dim3D = 3 ;
  GeometricModelPtr gm = 
    new GeometricModel(dim3D, region_list, (Epetra_MpiComm *)comm);
  
  // intantiate the mesh: build the geometric model
  string fname = string("./test/tmp_proto_fs.exo") ;
  Teuchos::RCP<Mesh> mesh = 
    Teuchos::rcp(new Mesh_MSTK(fname.c_str(), comm, gm));

  // create and initialize the deform mesh class
  DeformMesh deform_test(plist,mesh);
  deform_test.layer_profile();

  // say goodbye and exit
  deform_test.print_goodbye();
}
