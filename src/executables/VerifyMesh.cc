/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: William Perkins
*/

#include <fstream>

#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Kokkos_Core.hpp"

#include "AmanziComm.hh"
#include "MeshFactory.hh"
#include "MeshAudit.hh"
#include "MeshException.hh"
#include "VerboseObject_objs.hh"


int
main(int argc, char* argv[])
{
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  Kokkos::initialize();
  int status;
  {
    auto comm = Amanzi::getDefaultComm();
    const int nproc(comm->NumProc());
    const int me(comm->MyPID());

    // handle command line

    Teuchos::CommandLineProcessor CLP;

    CLP.setDocString("reads mesh file or file set and does a series of checks\n");

    std::vector<Amanzi::AmanziMesh::Framework> frameworks_avail;
    std::vector<const char*> frameworks_avail_names;

    if (Amanzi::AmanziMesh::framework_enabled(Amanzi::AmanziMesh::Framework::MSTK)) {
      frameworks_avail.push_back(Amanzi::AmanziMesh::Framework::MSTK);
      frameworks_avail_names.push_back(
        Amanzi::AmanziMesh::to_string(Amanzi::AmanziMesh::Framework::MSTK).c_str());
    }
    if (Amanzi::AmanziMesh::framework_enabled(Amanzi::AmanziMesh::Framework::MOAB)) {
      frameworks_avail.push_back(Amanzi::AmanziMesh::Framework::MOAB);
      frameworks_avail_names.push_back(
        Amanzi::AmanziMesh::to_string(Amanzi::AmanziMesh::Framework::MOAB).c_str());
    }
    if (frameworks_avail.size() == 0) {
      std::cerr << "error: this build has no frameworks which can read mesh files." << std::endl;
      return 1;
    }

    Amanzi::AmanziMesh::Framework the_framework = frameworks_avail[0];
    CLP.setOption("framework",
                  &the_framework,
                  frameworks_avail.size(),
                  frameworks_avail.data(),
                  frameworks_avail_names.data(),
                  "mesh framework preference",
                  false);

    std::string filename;
    CLP.setOption("file", &filename, "name of file or file set", true);

    bool dump_face_map(false);
    CLP.setOption("face-map", "no-face-map", &dump_face_map, "print the face Epetra_Map");

    bool dump_cell_map(false);
    CLP.setOption("cell-map", "no-cell-map", &dump_cell_map, "print the cell Epetra_Map");

    bool dump_node_map(false);
    CLP.setOption("node-map", "no-node-map", &dump_node_map, "print the node Epetra_Map");


    CLP.throwExceptions(false);

    int ierr(0), aerr(0);
    Teuchos::CommandLineProcessor::EParseCommandLineReturn parseReturn;
    try {
      parseReturn = CLP.parse(argc, argv);
    } catch (const std::exception& e) {
      std::cerr << "error: " << e.what() << std::endl;
      ierr++;
    }

    comm->SumAll(&ierr, &aerr, 1);

    if (aerr > 0) { return 1; }

    if (parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) { return 0; }

    if (parseReturn != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) { return 1; }

    // One command line argument is a file name. Three
    // types are supported depending on which frameworks are compiled in

    // A second command line argument is the framework preference.
    //
    // Put the user's preference first, but then allow all possible frameworks
    // in case they messed up.
    Amanzi::AmanziMesh::MeshFactory meshfactory(comm);
    std::vector<Amanzi::AmanziMesh::Framework> frameworks;
    frameworks.push_back(the_framework);
    for (auto p : frameworks_avail) {
      if (p != the_framework) frameworks.push_back(p);
    }
    Amanzi::AmanziMesh::Preference prefs(frameworks);


    Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh;

    ierr = 0;
    aerr = 0;
    try {
      if (me == 0) {
        std::cerr << "Attempting to read \"" << filename << "\" with ";
        for (auto p : frameworks) {
          std::cerr << "\"" << Amanzi::AmanziMesh::to_string(p) << "\",";
        }
        std::cerr << std::endl;
      }
      meshfactory.set_preference(prefs);
      mesh = meshfactory.create(filename);
    } catch (const Amanzi::AmanziMesh::Message& e) {
      std::cerr << argv[0] << ": mesh error: " << e.what() << std::endl;
      ierr++;
    } catch (const std::exception& e) {
      std::cerr << argv[0] << ": error: " << e.what() << std::endl;
      ierr++;
    }

    comm->SumAll(&ierr, &aerr, 1);

    if (aerr > 0) { return 3; }


    if (nproc == 1) {
      Amanzi::AmanziMesh::MeshAudit audit(mesh);
      status = audit.Verify();
    } else {
      std::ostringstream ofile;
      ofile << "mesh_audit_" << std::setfill('0') << std::setw(4) << me << ".txt";
      std::ofstream ofs(ofile.str().c_str());
      if (me == 0) std::cout << "Writing results to " << ofile.str() << ", etc." << std::endl;
      Amanzi::AmanziMesh::MeshAudit audit(mesh, ofs);
      status = audit.Verify();
    }

    if (me == 0) { std::cout << filename << ": " << (status ? "has errors" : "OK") << std::endl; }

    if (dump_node_map) {
      if (me == 0) { std::cout << "Node Epetra Map: " << std::endl; }
      (mesh->getMap(Amanzi::AmanziMesh::Entity_kind::NODE, false)).Print(std::cout);
    }

    if (dump_face_map) {
      if (me == 0) { std::cout << "Face Epetra Map: " << std::endl; }
      (mesh->getMap(Amanzi::AmanziMesh::Entity_kind::FACE, false)).Print(std::cout);
    }

    if (dump_cell_map) {
      if (me == 0) { std::cout << "Cell Epetra Map: " << std::endl; }
      (mesh->getMap(Amanzi::AmanziMesh::Entity_kind::CELL, false)).Print(std::cout);
    }
  }
  Kokkos::finalize();
  return status;
}
