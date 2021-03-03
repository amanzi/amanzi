#include <UnitTest++.h>
#include <TestReporterStdout.h>
#include <mpi.h>

#include <fstream>
#include "Epetra_Map.h"
#include "AmanziComm.hh"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_ParameterXMLFileReader.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_StandardParameterEntryValidators.hpp"
#include "Teuchos_VerboseObjectParameterListHelpers.hpp"

#include "MeshAudit.hh"

#include "../../mesh_factory/MeshFrameworkFactory.hh"
#include "../Mesh_MSTK.hh"


TEST(ELIM_DEGEN_INLINE_PARTITION)
{
  
  std::string xml_filename = "test/po_test_pri.xml";
  std::string out_exo_filename = "test/po_mesh_out.exo";
  
  auto comm = Amanzi::getDefaultComm();
  int num_procs = comm->NumProc();
  int rank = comm->MyPID();
  
  std::cout << "Reading the input file..." << std::endl;
  Teuchos::RCP<Teuchos::ParameterList> plist = Teuchos::getParametersFromXmlFile(xml_filename);
  // create the geometric model and regions
  Teuchos::ParameterList reg_params = plist->sublist("regions");
  std::cout << "Creating the geometric model..." << std::endl;
  Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel> gm =
  Teuchos::rcp(new Amanzi::AmanziGeometry::GeometricModel(3, reg_params, *comm));
  
  // create and register meshes
  Teuchos::ParameterList mesh_plist = plist->sublist("mesh");
  Amanzi::AmanziMesh::MeshFrameworkFactory meshfactory(comm, gm);
  Amanzi::AmanziMesh::Preference prefs(meshfactory.get_preference());
  prefs.clear();
  prefs.push_back(Amanzi::AmanziMesh::Framework::MSTK);
  
  // create the base mesh
  Teuchos::RCP<Amanzi::AmanziMesh::MeshFramework> mesh;
  
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
  mesh = meshfactory.create(in_exo_file);
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
      
      comm->SumAll(&ierr, &aerr, 1);
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
  
  // Create the surface mesh if needed
  Teuchos::RCP<Amanzi::AmanziMesh::MeshFramework> surface_mesh = Teuchos::null;
  Teuchos::RCP<Amanzi::AmanziMesh::MeshFramework> surface3D_mesh = Teuchos::null;
  if (mesh_plist.isSublist("surface mesh")) {
    Teuchos::ParameterList surface_plist = mesh_plist.sublist("surface mesh");
    
    std::vector<std::string> setnames;
    if (surface_plist.isParameter("surface sideset name")) {
      setnames.push_back(surface_plist.get<std::string>("surface sideset name"));
    } else if (surface_plist.isParameter("surface sideset names")) {
      setnames = surface_plist.get<Teuchos::Array<std::string> >("surface sideset names").toVector();
    } else {
      Errors::Message message("Surface mesh sublist missing parameter \"surface sideset names\".");
      Exceptions::amanzi_throw(message);
    }
    
    if (mesh->manifold_dimension() == 3) {
      surface3D_mesh = meshfactory.create(mesh,setnames,Amanzi::AmanziMesh::FACE,false);
      surface_mesh = meshfactory.create(mesh,setnames,Amanzi::AmanziMesh::FACE,true);
    } else {
      surface3D_mesh = mesh;
      surface_mesh = meshfactory.create(mesh,setnames,Amanzi::AmanziMesh::CELL,true);
    }
    
    bool surf_verify = surface_plist.get<bool>("verify mesh", false);
    if (surf_verify) {
      if (rank == 0)
        std::cout << "Verifying surface mesh with Mesh Audit..." << std::endl;
      if (num_procs == 1) {
        Amanzi::MeshAudit surf_mesh_auditor(surface_mesh);
        int status = surf_mesh_auditor.Verify();
        if (status == 0) {
          std::cout << "Mesh Audit confirms that surface mesh is ok" << std::endl;
        } else {
          Errors::Message msg("Mesh Audit could not verify correctness of surface mesh.");
          Exceptions::amanzi_throw(msg);
        }
      } else {
        std::ostringstream ofile;
        ofile << "surf_mesh_audit_" << std::setfill('0') << std::setw(4) << rank << ".txt";
        std::ofstream ofs(ofile.str().c_str());
        if (rank == 0)
          std::cout << "Writing Surface Mesh Audit output to " << ofile.str() << ", etc." << std::endl;
        
        int ierr = 0, aerr = 0;
        Amanzi::MeshAudit surf_mesh_auditor(mesh, ofs);
        int status = surf_mesh_auditor.Verify();        // check the mesh
        if (status != 0) ierr = 1;
        
        comm->SumAll(&ierr, &aerr, 1);
        if (aerr == 0) {
          if (rank == 0)
            std::cout << "Mesh Audit confirms that surface mesh is ok" << std::endl;
        } else {
          Errors::Message msg("Mesh Audit could not verify correctness of surface mesh.");
          Exceptions::amanzi_throw(msg);
        }
      }
      std::cout << "Verifying the surface mesh using the internal MSTK check..." << std::endl;
      Amanzi::AmanziMesh::Mesh_MSTK *surf_mstk_mesh =
        dynamic_cast<Amanzi::AmanziMesh::Mesh_MSTK *>(surface_mesh.get());
      CHECK(surf_mstk_mesh->run_internal_mstk_checks());
    }  // if surf_verify
    
    if (surface_plist.isParameter("export mesh to file")) {
      std::string export_surf_mesh_filename =
      surface_plist.get<std::string>("export mesh to file");
      surface3D_mesh->write_to_exodus_file(export_surf_mesh_filename);
    }
  }
  
  mesh->write_to_exodus_file(out_exo_filename);
}
