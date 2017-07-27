/* -*-  mode: c++; indent-tabs-mode: nil -*- */
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
#include "MeshSurfaceCell.hh"
#include "GeometricModel.hh"

#include "ats_mesh_factory.hh"

namespace ATS {

void
createMesh(Teuchos::ParameterList& mesh_plist,
           const Teuchos::RCP<Epetra_MpiComm>& comm,
           const Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel>& gm,
           Amanzi::State& S)
{
  auto mesh_type = mesh_plist.get<std::string>("mesh type");
  if (mesh_type == "read mesh file") {
    // create the MSTK factory
    Amanzi::AmanziMesh::MeshFactory factory(comm.get());
    Amanzi::AmanziMesh::FrameworkPreference prefs(factory.preference());
    prefs.clear();
    prefs.push_back(Amanzi::AmanziMesh::MSTK);
    
    // from file
    Teuchos::ParameterList read_params = mesh_plist.sublist("read mesh file parameters");

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
    auto mesh = factory.create(file, gm);
    bool deformable = mesh_plist.get<bool>("deformable mesh",false);

    checkVerifyMesh(mesh_plist, mesh);
    S.RegisterMesh(Amanzi::Keys::cleanPListName(mesh_plist.name()), mesh, deformable);
    
  } else if (mesh_type == "generate mesh") {
    // create the MSTK factory
    Amanzi::AmanziMesh::MeshFactory factory(comm.get());
    Amanzi::AmanziMesh::FrameworkPreference prefs(factory.preference());
    prefs.clear();
    prefs.push_back(Amanzi::AmanziMesh::MSTK);

    // generated mesh
    auto mesh = factory.create(mesh_plist.sublist("generate mesh parameters"), gm);
    bool deformable = mesh_plist.get<bool>("deformable mesh",false);

    checkVerifyMesh(mesh_plist, mesh);
    S.RegisterMesh(Amanzi::Keys::cleanPListName(mesh_plist.name()), mesh, deformable);

  } else if (mesh_type == "logical mesh") {
    // -- from logical mesh file
    Amanzi::AmanziMesh::MeshLogicalFactory fac(comm.get(), gm);
    auto mesh = fac.Create(mesh_plist.sublist("logical mesh parameters"));
    bool deformable = mesh_plist.get<bool>("deformable mesh",false);

    checkVerifyMesh(mesh_plist, mesh);
    S.RegisterMesh(Amanzi::Keys::cleanPListName(mesh_plist.name()), mesh, deformable);

  } else if (mesh_type == "aliased") {
    S.AliasMesh(mesh_plist.sublist("aliased parameters").get<std::string>("alias"), Amanzi::Keys::cleanPListName(mesh_plist.name()));
    
  } else if (mesh_type == "surface") {
    Teuchos::ParameterList& surface_plist = mesh_plist.sublist("surface parameters");
    std::vector<std::string> setnames;
    if (surface_plist.isParameter("surface sideset name")) {
      setnames.push_back(surface_plist.get<std::string>("surface sideset name"));
    } else if (surface_plist.isParameter("surface sideset names")) {
      setnames = surface_plist.get<Teuchos::Array<std::string> >("surface sideset names").toVector();
    } else {
      Errors::Message message("Surface mesh sublist missing parameter \"surface sideset names\".");
      Exceptions::amanzi_throw(message);
    }

    Teuchos::RCP<Amanzi::AmanziMesh::Mesh> surface3D_mesh = Teuchos::null;
    Teuchos::RCP<Amanzi::AmanziMesh::Mesh> surface_mesh = Teuchos::null;

    // create the MSTK factory
    Amanzi::AmanziMesh::MeshFactory factory(comm.get());
    Amanzi::AmanziMesh::FrameworkPreference prefs(factory.preference());
    prefs.clear();
    prefs.push_back(Amanzi::AmanziMesh::MSTK);

    auto parent = S.GetMesh(surface_plist.get<std::string>("parent domain", "domain"));
    if (parent->manifold_dimension() == 3) {
      surface3D_mesh = factory.create(&*parent,setnames,Amanzi::AmanziMesh::FACE,false,false);
      surface_mesh = factory.create(&*parent,setnames,Amanzi::AmanziMesh::FACE,true,false);
    } else {
      surface_mesh = factory.create(&*parent,setnames,Amanzi::AmanziMesh::CELL,true,false);
    }
    bool deformable = mesh_plist.get<bool>("deformable mesh",false);

    if (surface3D_mesh.get()) {
      S.RegisterMesh(Amanzi::Keys::cleanPListName(mesh_plist.name())+"_3d", surface3D_mesh, deformable);
    } else {
      S.AliasMesh(surface_plist.get<std::string>("parent domain", "domain"),
                  Amanzi::Keys::cleanPListName(mesh_plist.name())+"_3d");
    }
    checkVerifyMesh(mesh_plist, surface_mesh);
    S.RegisterMesh(Amanzi::Keys::cleanPListName(mesh_plist.name()), surface_mesh, deformable);

  } else if (mesh_type == "column") {
    Teuchos::ParameterList& column_list = mesh_plist.sublist("column parameters");
    int lid = column_list.get<int>("entity LID");
    auto parent = S.GetMesh(column_list.get<std::string>("parent domain", "domain"));
    auto mesh = Teuchos::rcp(new Amanzi::AmanziMesh::MeshColumn(*parent, lid));
    bool deformable = mesh_plist.get<bool>("deformable mesh",false);

    checkVerifyMesh(mesh_plist, mesh);
    S.RegisterMesh(Amanzi::Keys::cleanPListName(mesh_plist.name()), mesh, deformable);

  } else if (mesh_type == "column surface") {
    Teuchos::ParameterList& column_list = mesh_plist.sublist("column surface parameters");
    std::string surface_setname = column_list.get<std::string>("subgrid set name", "surface");
    auto parent = S.GetMesh(column_list.get<std::string>("parent domain"));
    auto mesh = Teuchos::rcp(new Amanzi::AmanziMesh::MeshSurfaceCell(*parent, surface_setname));

    bool deformable = mesh_plist.get<bool>("deformable mesh",false);

    checkVerifyMesh(mesh_plist, mesh);
    S.RegisterMesh(Amanzi::Keys::cleanPListName(mesh_plist.name()), mesh, deformable);

  } else if (mesh_type == "subgrid") {
    Teuchos::ParameterList& subgrid = mesh_plist.sublist("subgrid parameters");
    auto kind_str = subgrid.get<std::string>("entity kind");

    Amanzi::AmanziMesh::Entity_kind kind;
    if (kind_str == "CELL" || kind_str == "cell" || kind_str == "Cell") {
      kind = Amanzi::AmanziMesh::CELL;
    } else if (kind_str == "FACE" || kind_str == "face" || kind_str == "Face") {
      kind = Amanzi::AmanziMesh::FACE;
    } else if (kind_str == "NODE" || kind_str == "node" || kind_str == "Node") {
      kind = Amanzi::AmanziMesh::NODE;
    } else {
      Errors::Message msg("ATS Mesh Factory: Cannot create a subgrid mesh on the requested entity kind");
      Exceptions::amanzi_throw(msg);
    }

    auto region = subgrid.get<std::string>("region", "ENTIRE_MESH_REGION");
    auto parent_mesh = S.GetMesh(subgrid.get<std::string>("parent mesh", "domain"));
    bool flyweight = subgrid.get<bool>("flyweight mesh", false);
    
    auto comm_self = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_SELF));
    
    // for each id in the regions of the parent mesh on entity, create a subgrid mesh
    Amanzi::AmanziMesh::Entity_ID_List entities;
    parent_mesh->get_set_entities(region, kind, Amanzi::AmanziMesh::OWNED, &entities);
    const Epetra_Map& map = parent_mesh->map(kind,false);
    
    for (auto lid : entities) {
      Amanzi::AmanziMesh::Entity_ID gid = map.GID(lid);
      std::stringstream name;
      name << Amanzi::Keys::cleanPListName(mesh_plist.name()) << "_" << gid;

      Teuchos::ParameterList subgrid_i_list;
      if (subgrid.isSublist(name.str())) {
        subgrid_i_list = subgrid.sublist(name.str());
      } else {
        subgrid_i_list = subgrid.sublist(Amanzi::Keys::cleanPListName(mesh_plist.name())+"_*");
      }
      subgrid_i_list.setName(name.str());
      subgrid_i_list.set("entity kind", kind_str);
      subgrid_i_list.set("entity LID", lid);
      subgrid_i_list.set("subgrid set name", region);
      createMesh(subgrid_i_list, comm_self, gm, S);
    }

  } else {
    Errors::Message msg;
    msg << "ATS Mesh Factory: unknown \"mesh type\" parameter \"" << mesh_type
        << "\" in mesh \"" << Amanzi::Keys::cleanPListName(mesh_plist.name()) << "\".";
    Exceptions::amanzi_throw(msg);
  }
}


bool checkVerifyMesh(Teuchos::ParameterList& mesh_plist,
                     Teuchos::RCP<Amanzi::AmanziMesh::Mesh> mesh)
{
  // mesh verification
  ASSERT(!mesh.is_null());
  bool verify = mesh_plist.get<bool>("verify mesh", false);
  if (verify) {

    int num_procs = mesh->get_comm()->NumProc();
    int rank = mesh->get_comm()->MyPID();
    
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
        return false;
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
      
      mesh->get_comm()->SumAll(&ierr, &aerr, 1);
      if (aerr == 0) {
        if (mesh->get_comm()->MyPID() == 0)
          std::cout << "Mesh Audit confirms that mesh is ok" << std::endl;
      } else {
        Errors::Message msg("Mesh Audit could not verify correctness of mesh.");
        Exceptions::amanzi_throw(msg);
        return false;
      }
    }
  }  // if verify
  return true;
}


void
createMeshes(Teuchos::ParameterList& global_list,
             const Teuchos::RCP<Epetra_MpiComm>& comm,
             const Teuchos::RCP<Amanzi::AmanziGeometry::GeometricModel>& gm,
             Amanzi::State& S)
{
  Teuchos::RCP<Teuchos::Time> volmeshtime =
      Teuchos::TimeMonitor::getNewCounter("volume mesh creation");
  Teuchos::TimeMonitor timer(*volmeshtime);

  Teuchos::ParameterList& meshes_list = global_list.sublist("mesh");

  // always try to do the domain mesh first
  if (meshes_list.isSublist("domain")) {
    createMesh(meshes_list.sublist("domain"), comm, gm, S);
  }

  // always try to do the surface mesh second
  if (meshes_list.isSublist("surface")) {
    createMesh(meshes_list.sublist("surface"), comm, gm, S);
  }

  // now do the rest
  for (auto sublist : meshes_list) {
    if (sublist.first != "domain" && sublist.first != "surface" &&
        meshes_list.isSublist(sublist.first)) {
      createMesh(meshes_list.sublist(sublist.first), comm, gm, S);
    }
  }

  // FIXME --etc
  // this should be dealt with somewhere else, and more generally
  // generalize vis for columns
  if (global_list.isSublist("visualization columns")) {
    auto surface_mesh = S.GetMesh("surface");
    Teuchos::ParameterList& vis_ss_plist = global_list.sublist("visualization columns"); 
    int nc = surface_mesh->num_entities(Amanzi::AmanziMesh::CELL, Amanzi::AmanziMesh::OWNED);
    for (int c=0; c!=nc; ++c){
      int id = surface_mesh->cell_map(false).GID(c);
      std::stringstream name_ss;
      name_ss << "column_" << id;
      vis_ss_plist.set("file name base", "visdump_"+name_ss.str());         
      global_list.set("visualization " +name_ss.str(), vis_ss_plist);
    }  
    global_list.remove("visualization columns");
  }

  // generalize vis for surface columns
  if (global_list.isSublist("visualization surface cells")) {
    auto surface_mesh = S.GetMesh("surface");
    Teuchos::ParameterList& vis_sf_plist = global_list.sublist("visualization surface cells");
    int nc = surface_mesh->num_entities(Amanzi::AmanziMesh::CELL, Amanzi::AmanziMesh::OWNED);
    for (int c=0; c!=nc; ++c){
      int id = surface_mesh->cell_map(false).GID(c);
      std::stringstream name_ss, name_sf;
      name_sf << "column_" << id << "_surface";
      vis_sf_plist.set("file name base", "visdump_"+name_sf.str());
      global_list.set("visualization " +name_sf.str(), vis_sf_plist);
    }
    global_list.remove("visualization surface cells");
  }
  
  Teuchos::TimeMonitor::summarize();
  Teuchos::TimeMonitor::zeroOutTimers();
}


} // namespace ATS
