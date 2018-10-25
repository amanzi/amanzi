#include <UnitTest++.h>
#include <TestReporterStdout.h>
#include <mpi.h>

#include <fstream>
#include "Map_type.h"
#include "Teuchos_DefaultMpiComm.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"

#include "VerboseObject_objs.hh"

#include "MeshAudit.hh"

#include "../../mesh_factory/MeshFactory.hh"
#include "../Mesh_MSTK.hh"


TEST(ELIM_DEGEN_PREPARTITION)
{
  std::string xml_filename = "test/test_degen_prepart.xml";
  
  auto comm = Comm_ptr_type( new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
  int num_procs = comm_->getSize();
  int rank = comm_->getRank();
  
  std::cout << "Reading the input file..." << std::endl;
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xml_filename);
  // create the geometric model and regions
  Teuchos::ParameterList reg_params = plist->sublist("regions");
  std::cout << "Creating the geometric model..." << std::endl;
  Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm =
  Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, reg_params, comm_.get()));
  
  // create and register meshes
  Teuchos::ParameterList mesh_plist = plist->sublist("mesh");
  Amanzi::AmanziMesh::MeshFactory factory(comm_.get());
  Amanzi::AmanziMesh::FrameworkPreference prefs(factory.preference());
  prefs.clear();
  prefs.push_back(Amanzi::AmanziMesh::MSTK);
  
  // create the base mesh
  Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh;
  
  std::string in_exo_file;
  if (mesh_plist.isSublist("read mesh file")) {
    Teuchos::ParameterList read_params = mesh_plist.sublist("read mesh file");
    
    // file name
    if (read_params.isParameter("file")) {
      in_exo_file = read_params.get<std::string>("file");
    } else {
      Errors::Message msg("\"read mesh file\" list missing \"file\" parameter.");
      Exceptions::amanzi_throw(msg);
    }
    
    // file format
    std::string format;
    if (read_params.isParameter("format")) {
      // Is the format one that we can read?
      format = read_params.get<std::string>("format");
      if (format != "Exodus II") {
        Errors::Message msg;
        msg << "\"read mesh file\" parameter \"format\" with value \"" << format
        << "\" not understood: valid formats are: \"Exodus II\".";
        Exceptions::amanzi_throw(msg);
      }
    } else {
      Errors::Message msg("\"read mesh file\" parameter \"format\" missing.");
      Exceptions::amanzi_throw(msg);
    }
  } else {
    Errors::Message msg("Must specify mesh sublist of type: \"read mesh file\".");
    Exceptions::amanzi_throw(msg);
  }
  std::cout << "Reading the mesh..." << std::endl;
  mesh = factory.create(in_exo_file, gm);
  AMANZI_ASSERT(!mesh.is_null());
  
  // mesh verification
  bool verify = mesh_plist.get<bool>("verify mesh", false);
  if (verify) {
    if (rank == 0)
      std::cout << "Verifying mesh with Mesh Audit..." << std::endl;
    if (num_procs == 1) {
      Amanzi::MeshAudit mesh_auditor(mesh);
      int status = mesh_auditor.Verify();
      if (status == 0) {
        std::cout << "Mesh Audit confirms that mesh is ok" << std::endl;
      } else {
        Errors::Message msg("Mesh Audit could not verify correctness of mesh.");
        Exceptions::amanzi_throw(msg);
      }
    } else {
      std::ostringstream ofile;
      ofile << "mesh_audit_" << std::setfill('0') << std::setw(4) << rank << ".txt";
      std::ofstream ofs(ofile.str().c_str());
      if (rank == 0)
        std::cout << "Writing Mesh Audit output to " << ofile.str() << ", etc." << std::endl;
      
      int ierr = 0, aerr = 0;
      Amanzi::MeshAudit mesh_auditor(mesh, ofs);
      int status = mesh_auditor.Verify();        // check the mesh
      if (status != 0) ierr = 1;
      
      comm_->SumAll(&ierr, &aerr, 1);
      if (aerr == 0) {
        if (rank == 0)
          std::cout << "Mesh Audit confirms that mesh is ok" << std::endl;
      } else {
        Errors::Message msg("Mesh Audit could not verify correctness of mesh.");
        Exceptions::amanzi_throw(msg);
      }
    }
  }  // if verify
  
  std::cout << "Verifying the mesh using the internal MSTK check..." << std::endl;
  Amanzi::AmanziMesh::Mesh_MSTK *mstk_mesh = 
      dynamic_cast<Amanzi::AmanziMesh::Mesh_MSTK *>(mesh.get());
  CHECK(mstk_mesh->run_internal_mstk_checks());
}
