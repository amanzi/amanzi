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
  std::cout << "Test: DeformMesh using a simple quadrilateral mesh" << endl;

  // sequential/parallel communicator object
#ifdef HAVE_MPI
  Epetra_MpiComm *comm = new Epetra_MpiComm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm *comm = new Epetra_SerialComm();
#endif

  // read parameter list from input XML file
  Teuchos::ParameterList plist;
  string xml_fname = "test/simple_mesh_deform_test_2D.xml";
  updateParametersFromXmlFile(xml_fname, &plist);

  // intantiate the mesh: get the region list
  Teuchos::ParameterList region_list = plist.get<Teuchos::ParameterList>("Regions");
  
  // intantiate the mesh: build the geometric model
  const int dim2D = 2 ;
  GeometricModelPtr gm = new GeometricModel(dim2D, region_list, (Epetra_MpiComm *)comm);
  
  // intantiate the mesh: build the geometric model
  int nx(100), ny(100) ;
  Teuchos::RCP<Mesh> mesh = 
    Teuchos::rcp(new Mesh_MSTK( 0., 0., 1., 1., nx, ny, comm, gm));

  // create and initialize the deform mesh class
  DeformMesh deform_test(plist,mesh);

#if 0
  // VTK output, starting mesh
  deform_test.print_VTK_unstructured_mesh( string("mesh_0") );

  // check nodal coordinates
  deform_test.check_mesh_nodes();
  LINE(--);

  // move a single node (top, middle)
  //deform_test.move_a_single_node();
  deform_test.move_a_node_column();

  // check again the nodal coordinates
  deform_test.check_mesh_nodes();
  LINE(--);

  // VTK output, final mesh
  deform_test.print_VTK_unstructured_mesh( string("mesh_1") );
#endif

  deform_test.parabolic_profile();

  // say goodbye and exit
  deform_test.print_goodbye();
}
