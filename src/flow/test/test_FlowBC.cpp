#include <iostream>
#include "UnitTest++.h"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

#include "Epetra_MpiComm.h"
#include "Epetra_SerialComm.h"

#include "Mesh_simple.hh"
#include "FlowBC.hpp"

TEST(FlowBC) {

#ifdef HAVE_MPI
  Epetra_MpiComm *comm = new Epetra_MpiComm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm *comm = new Epetra_SerialComm();
#endif

  // Create a simple 2x2x2 brick mesh
  Teuchos::RCP<Mesh> mesh(new Mesh_simple(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2, 2, 2, comm));

  // Read the flow BC parameters from an XML file.
  Teuchos::ParameterList flow_BCs_list;
  std::string xmlInFileName = "test/test_FlowBC.xml";
  Teuchos::updateParametersFromXmlFile(xmlInFileName, &flow_BCs_list);

  // Create the flow BC object using the read parameters.
  FlowBC bc(flow_BCs_list, mesh);

  // Write out what we've got.
  for (int i = 0; i < bc.NumBC(); ++i) {
    std::cout << "Type...     " << bc[i].Type << std::endl;
    std::cout << "SetID...    " << bc[i].SetID << std::endl;
    std::cout << "Faces...    ";
    for (int j = 0; j < bc[i].Faces.size(); ++j) std::cout << " " << bc[i].Faces[j];
    std::cout << std::endl;
    if (bc[i].Type == FlowBC::PRESSURE_CONSTANT)
      std::cout << "Value...    " << bc[i].Value << std::endl;
    std::cout << std::endl;
  }
  
  CHECK_EQUAL(bc.NumBC(),6);
}


