//! OutputXDMF: writes XDMF+H5 files for VisIt
/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon
*/

/*
  XDMF implementation of an Output object.
*/

#include "OutputUtils.hh"
#include "OutputXDMF.hh"

namespace Amanzi {

OutputXDMF::OutputXDMF(Teuchos::ParameterList& plist,
		       const Teuchos::RCP<const AmanziMesh::Mesh>& mesh)
    : filenamebase_(plist.get<std::string>("file name base", "amanzi_vis")),
      mesh_(mesh),
      is_dynamic_(plist.get<bool>("dynamic mesh", false)),
      init_(false)
{}


// open and close files
void
OutputXDMF::InitializeCycle(double time, int cycle) {
  time_ = time;
  cycle_ = cycle;

  if (!init_) {
    // create and set up the FileHDF5 object
    h5_data_ = Teuchos::rcp(new FileHDF5(mesh_->get_comm(), filenamebase_+"_data.h5", FILE_CREATE));
  
    // create and set up the FileHDF5Mesh h5
    h5_mesh_ = Teuchos::rcp(new FileHDF5(mesh_->get_comm(), filenamebase_+"_mesh.h5", FILE_CREATE));

    // write the mesh
    auto global_sizes = WriteMesh_(cycle);
    write_mesh_ = true;

    // create and set up the FileXDMF object
    xdmf_ = Teuchos::rcp(new FileXDMF(filenamebase_, std::get<0>(global_sizes), std::get<1>(global_sizes), std::get<2>(global_sizes)));
    xdmf_->CreateTimestep(time, cycle, true);

    init_ = true;

  } else if (is_dynamic_) {
    WriteMesh_(cycle);
    write_mesh_ = true;

    xdmf_->CreateTimestep(time, cycle, true);
    h5_data_->OpenFile();

  } else {
    write_mesh_ = false;

    xdmf_->CreateTimestep(time, cycle, false);
    h5_data_->OpenFile();
  }    
  
  //h5_data_->CreateGroup(cycle_group_);
}


void OutputXDMF::FinalizeCycle() {
  h5_data_->CloseFile();
  xdmf_->CloseTimestep(time_, cycle_, write_mesh_);
}


std::string OutputXDMF::Filename() const
{
  Errors::Message message("OutputXDMF: This is a multi-part file, no single filename exists.");
  throw(message);
}



std::tuple<int,int,int>
OutputXDMF::WriteMesh_(int cycle)
{
  if (mesh_->is_logical()) {
    Errors::Message message("OutputXDMF: Mesh_Logical cannot yet be written (try ??? instead).");
    throw(message);
  }

  std::stringstream h5path;
  h5path << "/" << cycle << "/Mesh/";

  h5_mesh_->OpenFile();
  h5_mesh_->CreateGroup(h5path.str());
  
  // get node and cell maps
  auto node_map = mesh_->map(AmanziMesh::Entity_kind::NODE, false);
  auto global_nodes_tp = node_map->getGlobalNumElements();

  auto cell_map = mesh_->map(AmanziMesh::Entity_kind::CELL, false);
  Tpetra::global_size_t global_cells_tp = cell_map->getGlobalNumElements();
  
  int space_dim = mesh_->space_dimension();
  int manifold_dim = mesh_->manifold_dimension();
  
  // count total connections
  Tpetra::global_size_t local_conns(0); // length of MixedElements
  for (AmanziMesh::Entity_ID c=0; c!=cell_map->getNodeNumElements(); ++c) {
    AmanziMesh::Cell_type ctype = mesh_->cell_get_type(c);
    if (ctype != AmanziMesh::POLYGON) {
      Kokkos::View<AmanziMesh::Entity_ID*> nodes;
      mesh_->cell_get_nodes(c,nodes);
      local_conns += nodes.extent(0) + 1;

    } else if (manifold_dim == 2) {
      Kokkos::View<AmanziMesh::Entity_ID*> nodes;
      mesh_->cell_get_nodes(c,nodes);
      local_conns += nodes.extent(0) + 2;

    } else {
      Errors::Message message("OutputXDMF: Polyhedral meshes cannot yet be visualized in XDMF (try SILO instead).");
      throw(message);
    }
  }

  Tpetra::global_size_t global_conns_tp;
  Teuchos::reduceAll(*mesh_->get_comm(), Teuchos::REDUCE_SUM, 1, &local_conns, &global_conns_tp);

  // check that sizes aren't too big!  ASCEMIO can only handle things representable as INT
  // convert from Tpetra::global_size_t to int, failing if too big.
  if (global_nodes_tp > std::numeric_limits<int>::max() ||
      global_cells_tp > std::numeric_limits<int>::max() ||
      global_conns_tp > std::numeric_limits<int>::max()) {
    Errors::Message message("OutputXDMF: ASCEMIO does not support meshes with global sizes larger than a (signed) int.");
    throw(message);
  }
  int global_nodes = static_cast<int>(global_nodes_tp);
  int global_cells = static_cast<int>(global_cells_tp);
  int global_conns = static_cast<int>(global_conns_tp);
  
  // Get and write coordinates and coordinate map
  // -- create and store a vector of coordinates
  Kokkos::View<double**, Kokkos::LayoutRight, AmanziDefaultHost> coords("coordinates", node_map->getNodeNumElements(), 3);
  AmanziGeometry::Point p;
  for (AmanziMesh::Entity_ID i=0; i!=coords.extent(0); ++i) {
    mesh_->node_get_coordinates(i, &p);
    for (int j=0; j!=p.dim(); ++j) {
      coords(i,j) = p[j];
    }
  }
  // -- write the coordinates
  h5_mesh_->WriteView(coords, global_nodes, h5path.str()+"Nodes");

  // -- write the coordinate map
  h5_mesh_->WriteVector(*node_map, h5path.str()+"NodeMap");

  // Get and write connectivity information
  // nodes are written to h5 out of order, need the natural node map
  auto ghosted_natural_nodes = GetNaturalMap(mesh_->map(AmanziMesh::Entity_kind::NODE, true),
          mesh_->map(AmanziMesh::Entity_kind::NODE, false));
  
  // Create a map and vector to store connection info
  auto map = Teuchos::rcp(new Map_type(global_conns, local_conns, 0, mesh_->get_comm()));
  IntVector_type conns(map);
  {
    int lcv = 0;
    int lcv_entity = 0;
    auto connv = conns.getLocalViewHost();
    
    for (int c=0; c!=cell_map->getNodeNumElements(); ++c) {
      AmanziMesh::Cell_type ctype = mesh_->cell_get_type(c);
      if (ctype != AmanziMesh::POLYGON) {
        // store cell type id
        connv(lcv++,0) = XDMFCellTypeID(ctype);

        // store nodes in the correct order
        Kokkos::View<AmanziMesh::Entity_ID*> nodes;
        mesh_->cell_get_nodes(c, nodes);

        for (int i=0; i!=nodes.extent(0); ++i) {
          connv(lcv++, 0) = ghosted_natural_nodes->getGlobalElement(nodes(i));
        }

      } else if (space_dim == 2) {
        // store cell type id
        connv(lcv++,0) = XDMFCellTypeID(ctype);

        // store node count, then nodes in the correct order
        Kokkos::View<AmanziMesh::Entity_ID*> nodes;
        mesh_->cell_get_nodes(c, nodes);
        connv(lcv++,0) = nodes.extent(0);

        for (int i=0; i!=nodes.extent(0); ++i) {
          connv(lcv++,0) = ghosted_natural_nodes->getGlobalElement(nodes(i));
        }
      }
    }
  }

  // write the connections
  h5_mesh_->WriteVector(conns, h5path.str()+"MixedElements");

  // Write the cell map
  h5_mesh_->WriteVector(*cell_map, h5path.str()+"ElementMap");
  
  h5_mesh_->CloseFile();
  return std::make_tuple(global_nodes, global_cells, global_conns);  
}


} // namespace Amanzi
