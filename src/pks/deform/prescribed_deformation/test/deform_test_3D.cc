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
  Teuchos::Ptr<Teuchos::ParameterList> plist;
  string xml_fname = "test/simple_mesh_deform_test_3D.xml";
  updateParametersFromXmlFile(xml_fname, plist);

  // intantiate the mesh: get the region list
  Teuchos::ParameterList& region_list = 
    plist->get<Teuchos::ParameterList>("Regions");
  
  // intantiate the mesh: build the geometric model
  const int dim3D = 3 ;
  GeometricModelPtr gm = 
    new GeometricModel(dim3D, region_list, (Epetra_MpiComm *)comm);
  
  // intantiate the mesh: build the geometric model (big domain by default)
  string fname = string("./test/tmp_proto_fs.exo") ;
  
  // select a small subset of the domain
  bool do_subset_flag = false ;
  if ( do_subset_flag ) {
    fname = string("./test/poly0_subset_proto_fs.exo");
  }
  
  // this version uses three meshes:
  // mesh0 : starting mesh
  // mesh1 : final mesh (which is input from file fname)
  // mesh  : current mesh
  Teuchos::RCP<Mesh> mesh0 = 
    Teuchos::rcp(new Mesh_MSTK(fname.c_str(), comm, gm));

  Teuchos::RCP<Mesh> mesh1 = 
    Teuchos::rcp(new Mesh_MSTK(fname.c_str(), comm, gm));

  Teuchos::RCP<Mesh> mesh  = 
    Teuchos::rcp(new Mesh_MSTK(fname.c_str(), comm, gm));

  // create and initialize the deform mesh class
  DeformMesh deform_test(*plist,mesh0,mesh1,mesh);
  
  // get the list of nodes on the top
  //deform_test.mesh_deformation();
  //deform_test.analyze_final_mesh();

  // say goodbye and exit
  deform_test.print_goodbye();
}
