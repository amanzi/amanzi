/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
//! Simple wrapper that takes a ParameterList and generates all needed meshes.

/*
  ATS is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "Epetra_MpiComm.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_TimeMonitor.hpp"

#include "MeshAudit.hh"
#include "MeshFactory.hh"
#include "MeshLogicalFactory.hh"
#include "MeshColumn.hh"
#include "GeometricModel.hh"

#include "ats_mesh_factory.hh"

namespace ATS {

void
createMeshes(Teuchos::ParameterList& mesh_plist,
             const Teuchos::RCP<Epetra_MpiComm>& comm,
             const Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel>& gm,
             Amanzi::State& S) {

  // create the MSTK factory
  Amanzi::AmanziMesh::MeshFactory factory(comm.get());
  Amanzi::AmanziMesh::FrameworkPreference prefs(factory.preference());
  prefs.clear();
  prefs.push_back(Amanzi::AmanziMesh::MSTK);

  // create the base mesh
  Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh;

  if (mesh_plist.isSublist("read mesh file")) {
    // from file
    Teuchos::ParameterList read_params = mesh_plist.sublist("read mesh file");

    // file name
    std::string file;
    if (read_params.isParameter("file")) {
      file = read_params.get<std::string>("file");
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

    // create the mesh from the file
    Teuchos::RCP<Teuchos::Time> volmeshtime = Teuchos::TimeMonitor::getNewCounter("volume mesh creation");
    Teuchos::TimeMonitor timer(*volmeshtime);
    mesh = factory.create(file, gm);


  } else if (mesh_plist.isSublist("generate mesh")) {
    // generated mesh
    Teuchos::ParameterList gen_params = mesh_plist.sublist("generate mesh");
    mesh = factory.create(gen_params, gm);

  } else if (mesh_plist.isSublist("logical mesh")) {
    // -- from logical mesh file
    Amanzi::AmanziMesh::MeshLogicalFactory fac(comm.get(), gm);
    mesh = fac.Create(mesh_plist.sublist("logical mesh"));
    
  } else if (mesh_plist.isSublist("embedded logical mesh")) {
    Errors::Message msg("\"embedded logical mesh\" option not yet implemented.");
    Exceptions::amanzi_throw(msg);
    
  } else {
    Errors::Message msg("Must specify mesh sublist of type: \"read mesh file\", \"generate mesh\", or \"logical mesh\".");
    Exceptions::amanzi_throw(msg);
  }

  ASSERT(!mesh.is_null());

  // mesh verification
  bool verify = mesh_params.get<bool>("verify mesh", false);
  if (verify) {
    std::cerr << "Verifying mesh with Mesh Audit..." << std::endl;
    if (size == 1) {
      Amanzi::MeshAudit mesh_auditor(mesh);
      int status = mesh_auditor.Verify();
      if (status == 0) {
        std::cout << "Mesh Audit confirms that mesh is ok" << std::endl;
      } else {
        std::cout << "Mesh Audit could not verify correctness of mesh" << std::endl;
        return 1;
      }
    } else {
      std::ostringstream ofile;
      ofile << "mesh_audit_" << std::setfill('0') << std::setw(4) << rank << ".txt";
      std::ofstream ofs(ofile.str().c_str());
      if (rank == 0)
        std::cerr << "Writing Mesh Audit output to " << ofile.str() << ", etc." << std::endl;
      
      ierr = 0;
      Amanzi::MeshAudit mesh_auditor(mesh, ofs);
      int status = mesh_auditor.Verify();        // check the mesh
      if (status != 0) ierr++;
      
      comm->SumAll(&ierr, &aerr, 1);
      if (aerr == 0) {
        if (rank == 0)
          std::cerr << "Mesh Audit confirms that mesh is ok" << std::endl;
      } else {
        if (rank == 0)
          std::cerr << "Mesh Audit could not verify correctness of mesh" << std::endl;
        return 1;
      }
    }
  }  // if verify


  // Create the surface mesh if needed
  Teuchos::RCP<Amanzi::AmanziMesh::Mesh> surface_mesh = Teuchos::null;
  Teuchos::RCP<Amanzi::AmanziMesh::Mesh> surface3D_mesh = Teuchos::null;
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
      surface3D_mesh = factory.create(&*mesh,setnames,Amanzi::AmanziMesh::FACE,false,false);
      surface_mesh = factory.create(&*mesh,setnames,Amanzi::AmanziMesh::FACE,true,false);
    } else {
      surface3D_mesh = mesh;
      surface_mesh = factory.create(&*mesh,setnames,Amanzi::AmanziMesh::CELL,true,false);
    }

    bool surf_verify = surface_plist.get<bool>("verify mesh", false);
    if (surf_verify) {
      std::cerr << "Verifying surface mesh with Mesh Audit..." << std::endl;
      if (size == 1) {
        Amanzi::MeshAudit surf_mesh_auditor(surface_mesh);
        int status = surf_mesh_auditor.Verify();
        if (status == 0) {
          std::cout << "Mesh Audit confirms that surface mesh is ok" << std::endl;
        } else {
          std::cout << "Mesh Audit could not verify correctness of surface mesh" << std::endl;
          return 1;
        }
      } else {
        std::ostringstream ofile;
        ofile << "surf_mesh_audit_" << std::setfill('0') << std::setw(4) << rank << ".txt";
        std::ofstream ofs(ofile.str().c_str());
        if (rank == 0)
          std::cerr << "Writing Surface Mesh Audit output to " << ofile.str() << ", etc." << std::endl;
        
        ierr = 0;
        Amanzi::MeshAudit surf_mesh_auditor(mesh, ofs);
        int status = surf_mesh_auditor.Verify();        // check the mesh
        if (status != 0) ierr++;
        
        comm->SumAll(&ierr, &aerr, 1);
        if (aerr == 0) {
          std::cerr << "Surface Mesh Audit confirms that mesh is ok" << std::endl;
        } else {
          if (rank == 0)
            std::cerr << "Surface Mesh Audit could not verify correctness of mesh" << std::endl;
          return 1;
        }
      }
    }  // if surf_verify


    if (surface_plist.isParameter("export mesh to file")) {
      std::string export_surf_mesh_filename =
            surface_plist.get<std::string>("export mesh to file");
      surface3D_mesh->write_to_exodus_file(export_surf_mesh_filename);
    }
  }

  // column meshes
  std::vector<Teuchos::RCP<Amanzi::AmanziMesh::Mesh> > col_meshes;
  if (mesh_plist.isSublist("column meshes")) {
    int nc = mesh->num_columns();
    col_meshes.resize(nc, Teuchos::null);
    for (int c=0; c!=nc; ++c) {
      col_meshes[c] = Teuchos::rcp(new Amanzi::AmanziMesh::MeshColumn(*mesh, c));
    }
  }  
  
  Teuchos::TimeMonitor::summarize();
  Teuchos::TimeMonitor::zeroOutTimers();

  // register meshes with state
  bool deformable = mesh_plist.get<bool>("deformable mesh",false);
  S.RegisterDomainMesh(mesh, deformable);
  if (surface3D_mesh != Teuchos::null)
    S.RegisterMesh("surface_3d", surface3D_mesh, deformable);
  if (surface_mesh != Teuchos::null)
    S.RegisterMesh("surface", surface_mesh, deformable);

  if (col_meshes.size() > 0) {
    for (int c=0; c!=col_meshes.size(); ++c) {
      std::stringstream namestream;
      namestream << "column_" << c;
      S.RegisterMesh(namestream.str(), col_meshes[c], deformable);
    }
  }
}


} // namespace ATS
