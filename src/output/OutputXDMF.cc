/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

//! OutputXDMF: writes XDMF+H5 files for visualization in VisIt
/*
  XDMF implementation of an Output object, can only work as a Vis object as it
  needs a mesh and cannot handle face DoFs.
*/

#include "OutputUtils.hh"
#include "OutputXDMF.hh"

namespace Amanzi {

OutputXDMF::OutputXDMF(Teuchos::ParameterList& plist,
                       const Teuchos::RCP<const AmanziMesh::Mesh>& mesh)
  : mesh_(mesh), is_dynamic_(plist.get<bool>("dynamic mesh", false)), init_(false)
{
  filenamebase_ = plist.get<std::string>("file name base");
}


// open and close files
void
OutputXDMF::createTimestep(double time, int cycle)
{
  time_ = time;
  cycle_ = cycle;

  if (!init_) {
    // create and set up the FileHDF5 object
    h5_data_filename_ = filenamebase_ + "_data.h5";
    h5_data_ = std::make_unique<FileHDF5>(mesh_->getComm(), h5_data_filename_, FILE_CREATE);

    // create and set up the FileHDF5Mesh h5
    h5_mesh_filename_ = filenamebase_ + "_mesh.h5";
    h5_mesh_ = std::make_unique<FileHDF5>(mesh_->getComm(), h5_mesh_filename_, FILE_CREATE);

    // write the mesh
    auto global_sizes = writeMesh_(cycle);
    write_mesh_ = true;

    // create and set up the FileXDMF object
    xdmf_filename_ = filenamebase_;
    xdmf_ = std::make_unique<FileXDMF>(xdmf_filename_,
                                       std::get<0>(global_sizes),
                                       std::get<1>(global_sizes),
                                       std::get<2>(global_sizes));
    xdmf_->createTimestep(time, cycle, true);

    init_ = true;

  } else if (is_dynamic_) {
    writeMesh_(cycle);
    write_mesh_ = true;

    xdmf_->createTimestep(time, cycle, true);
    h5_data_->openFile();

  } else {
    write_mesh_ = false;

    xdmf_->createTimestep(time, cycle, false);
    h5_data_->openFile();
  }
  // h5_data_->createGroup(cycle_group_);
}


void
OutputXDMF::finalizeTimestep()
{
  h5_data_->closeFile();
  xdmf_->finalizeTimestep(time_, cycle_, write_mesh_);
}


std::tuple<int, int, int>
OutputXDMF::writeMesh_(int cycle)
{
  if (mesh_->isLogical()) return writeDualMesh_(cycle);
  auto vis_mesh = mesh_->getVisMesh();

  std::stringstream h5path;
  h5path << "/" << cycle << "/Mesh/";

  h5_mesh_->openFile();
  h5_mesh_->createGroup(h5path.str());

  // get node and cell maps
  auto node_map = vis_mesh.getMap(AmanziMesh::Entity_kind::NODE, false);
  auto global_nodes_tp = node_map->getGlobalNumElements();

  auto cell_map = vis_mesh.getMap(AmanziMesh::Entity_kind::CELL, false);
  Tpetra::global_size_t global_cells_tp = cell_map->getGlobalNumElements();

  int space_dim = vis_mesh.getSpaceDimension();
  int manifold_dim = vis_mesh.getManifoldDimension();

  // count total connections
  Tpetra::global_size_t local_conns(0); // length of MixedElements
  for (AmanziMesh::Entity_ID c = 0; c != cell_map->getLocalNumElements(); ++c) {
    AmanziMesh::Cell_kind ctype = vis_mesh.getCellKind(c);
    if (ctype == AmanziMesh::Cell_kind::POLYGON) {
      int nnodes = vis_mesh.getCellNodes(c).size();
      local_conns += nnodes + 2;

    } else if (ctype == AmanziMesh::Cell_kind::POLYHED) {
      int my_conns = 2; // POLYHED indicator + nfaces
      auto faces = vis_mesh.getCellFaces(c);
      for (auto f : faces) {
        my_conns += vis_mesh.getFaceNodes(f).size() + 1;
      }
      local_conns += my_conns;

    } else {
      int nnodes = vis_mesh.getCellNodes(c).size();
      local_conns += nnodes + 1;
    }
  }

  Tpetra::global_size_t global_conns_tp;
  Teuchos::reduceAll(*vis_mesh.getComm(), Teuchos::REDUCE_SUM, 1, &local_conns, &global_conns_tp);

  // check that sizes aren't too big!  ASCEMIO can only handle things
  // representable as INT convert from Tpetra::global_size_t to int, failing if
  // too big.
  if (global_nodes_tp > std::numeric_limits<int>::max() ||
      global_cells_tp > std::numeric_limits<int>::max() ||
      global_conns_tp > std::numeric_limits<int>::max()) {
    Errors::Message message("OutputXDMF: ASCEMIO does not support meshes with "
                            "global sizes larger than a (signed) int.");
    throw(message);
  }
  int global_nodes = static_cast<int>(global_nodes_tp);
  int global_cells = static_cast<int>(global_cells_tp);
  int global_conns = static_cast<int>(global_conns_tp);

  // Get and write coordinates and coordinate map
  // -- create and store a vector of coordinates
  Kokkos::View<double**, Kokkos::LayoutRight, DefaultHost> coords(
    "coordinates", node_map->getLocalNumElements(), 3);
  AmanziGeometry::Point p;
  for (AmanziMesh::Entity_ID i = 0; i != coords.extent(0); ++i) {
    p = vis_mesh.getNodeCoordinate(i);
    for (int j = 0; j != p.dim(); ++j) { coords(i, j) = p[j]; }
  }
  // -- write the coordinates
  h5_mesh_->writeView(h5path.str() + "Nodes", coords);

  // -- write the coordinate map
  h5_mesh_->writeVector(h5path.str() + "NodeMap", OutputUtils::asVector(node_map));

  // Get and write connectivity information
  // nodes are written to h5 out of order, need the natural node map
  auto [unused, ghosted_natural_nodes] =
    AmanziMesh::createMapsFromContiguousGIDs(vis_mesh, AmanziMesh::Entity_kind::NODE);

  // create a map and vector to store connection info
  auto map = Teuchos::rcp(new Map_type(global_conns, local_conns, 0, vis_mesh.getComm()));
  Vector_type_<GO> conns(map);
  {
    int lcv = 0;
    int lcv_entity = 0;
    auto connv = conns.getLocalViewHost(Tpetra::Access::ReadWrite);

    for (int c = 0; c != cell_map->getLocalNumElements(); ++c) {
      AmanziMesh::Cell_kind ctype = vis_mesh.getCellKind(c);
      if (ctype == AmanziMesh::Cell_kind::POLYHED) {
        std::cout << "Writing conn of type POLYHEDRON" << std::endl;
        auto faces = vis_mesh.getCellFaces(c);
        connv(lcv++, 0) = XDMFCellTypeID(ctype);
        connv(lcv++, 0) = faces.size();

        for (auto f : faces) {
          // store node count, then nodes in the correct order
          auto nodes = vis_mesh.getFaceNodes(f);
          connv(lcv++, 0) = nodes.size();

          for (int i = 0; i != nodes.size(); ++i) {
            connv(lcv++, 0) = ghosted_natural_nodes->getGlobalElement(nodes(i));
          }
        }

      } else if (ctype == AmanziMesh::Cell_kind::POLYGON) {
        std::cout << "Writing conn of type POLYGON" << std::endl;
        // store cell type id
        connv(lcv++, 0) = XDMFCellTypeID(ctype);

        // store node count, then nodes in the correct order
        auto nodes = vis_mesh.getCellNodes(c);
        connv(lcv++, 0) = nodes.size();

        for (int i = 0; i != nodes.size(); ++i) {
          connv(lcv++, 0) = ghosted_natural_nodes->getGlobalElement(nodes[i]);
        }

      } else {
        // store cell type id
        connv(lcv++, 0) = XDMFCellTypeID(ctype);

        // store nodes in the correct order
        auto nodes = vis_mesh.getCellNodes(c);

        for (int i = 0; i != nodes.size(); ++i) {
          connv(lcv++, 0) = ghosted_natural_nodes->getGlobalElement(nodes[i]);
        }


    }
    }
  }

  // write the connections
  h5_mesh_->writeVector(h5path.str() + "MixedElements", conns);

  // write the cell map
  h5_mesh_->writeVector(h5path.str() + "ElementMap", OutputUtils::asVector(cell_map));

  h5_mesh_->closeFile();
  return std::make_tuple(global_nodes, global_cells, global_conns);
}


std::tuple<int, int, int>
OutputXDMF::writeDualMesh_(int cycle)
{
  const auto& vis_mesh = mesh_->getVisMesh();

  std::stringstream h5path;
  h5path << "/" << cycle << "/Mesh/";

  h5_mesh_->openFile();
  h5_mesh_->createGroup(h5path.str());

  // on the dual, the node map is CELLS
  auto node_map = vis_mesh.getMap(AmanziMesh::Entity_kind::CELL, false);
  auto global_nodes_tp = node_map->getGlobalNumElements();

  int space_dim = vis_mesh.getSpaceDimension();
  int manifold_dim = vis_mesh.getManifoldDimension();

  // For the dual, connections are all faces with > 1 cells.
  // count total connections
  std::array<GO, 2> local_conn_ents = { 0, 0 };
  auto face_map = vis_mesh.getMap(AmanziMesh::Entity_kind::FACE, false);
  for (AmanziMesh::Entity_ID f = 0; f != face_map->getLocalNumElements(); ++f) {
    auto nfcells = vis_mesh.getFaceNumCells(f);
    if (nfcells > 1) {
      local_conn_ents[0] += nfcells;
      local_conn_ents[1]++;
    }
  }

  std::array<GO, 2> global_conn_ents;
  Teuchos::reduceAll(
    *vis_mesh.getComm(), Teuchos::REDUCE_SUM, 2, local_conn_ents.data(), global_conn_ents.data());
  GO global_conns_tp = global_conn_ents[0];
  GO global_cells_tp = global_conn_ents[1];

  // check that sizes aren't too big!  ASCEMIO can only handle things
  // representable as INT convert from Tpetra::global_size_t to int, failing if
  // too big.
  if (global_nodes_tp > std::numeric_limits<int>::max() ||
      global_cells_tp > std::numeric_limits<int>::max() ||
      global_conns_tp > std::numeric_limits<int>::max()) {
    Errors::Message message("OutputXDMF: ASCEMIO does not support meshes with "
                            "global sizes larger than a (signed) int.");
    throw(message);
  }
  int global_nodes = static_cast<int>(global_nodes_tp);
  int global_cells = static_cast<int>(global_cells_tp);
  int global_conns = static_cast<int>(global_conns_tp);

  // Get and write coordinates and coordinate map
  // -- create and store a vector of coordinates
  Kokkos::View<double**, Kokkos::LayoutRight, DefaultHost> coords(
    "coordinates", node_map->getLocalNumElements(), 3);
  AmanziGeometry::Point p;
  for (AmanziMesh::Entity_ID i = 0; i != coords.extent(0); ++i) {
    p = vis_mesh.getNodeCoordinate(i);
    for (int j = 0; j != p.dim(); ++j) { coords(i, j) = p[j]; }
  }
  // -- write the coordinates
  h5_mesh_->writeView(h5path.str() + "Nodes", coords);

  // -- write the coordinate map
  h5_mesh_->writeVector(h5path.str() + "NodeMap", OutputUtils::asVector(node_map));

  // Get and write connectivity information
  // nodes are written to h5 out of order, need the natural node map
  auto [unused, ghosted_natural_nodes] =
    AmanziMesh::createMapsFromContiguousGIDs(vis_mesh, AmanziMesh::Entity_kind::NODE);

  // create a map and vector to store connection info
  auto map = Teuchos::rcp(new Map_type(global_conns, local_conn_ents[0], 0, vis_mesh.getComm()));
  Vector_type_<GO> conns(map);
  {
    int lcv = 0;
    int lcv_entity = 0;
    auto connv = conns.getLocalViewHost(Tpetra::Access::ReadWrite);

    for (int f = 0; f != face_map->getLocalNumElements(); ++f) {
      auto cells = vis_mesh.getFaceCells(f);
      if (cells.size() > 1) {
        for (int i = 0; i != cells.size(); ++i) {
          connv(lcv++, 0) = ghosted_natural_nodes->getGlobalElement(cells[i]);
        }
      }
    }
  }

  // write the connections
  h5_mesh_->writeVector(h5path.str() + "MixedElements", conns);

  // write the cell map
  auto dual_cell_map =
    Teuchos::rcp(new Map_type(global_cells_tp, local_conn_ents[1], 0, vis_mesh.getComm()));
  h5_mesh_->writeVector(h5path.str() + "ElementMap", OutputUtils::asVector(dual_cell_map));
  h5_mesh_->closeFile();
  return std::make_tuple(global_nodes, global_cells, global_conns);
}


void
OutputXDMF::write(const Teuchos::ParameterList& attrs, const Vector_type& vec) const
{
  auto location = attrs.get<AmanziMesh::Entity_kind>("location");
  if (location == AmanziMesh::Entity_kind::CELL || location == AmanziMesh::Entity_kind::NODE) {
    xdmf_->writeField<Vector_type::scalar_type>(attrs.name(), location);
    std::stringstream path;
    path << "/" << attrs.name() << "/" << cycle_;
    h5_data_->writeVector<Vector_type::scalar_type>(path.str(), vec);
  }
}

void
OutputXDMF::write(const Teuchos::ParameterList& attrs, const IntVector_type& vec) const
{
  auto location = attrs.get<AmanziMesh::Entity_kind>("location");
  if (location == AmanziMesh::Entity_kind::CELL || location == AmanziMesh::Entity_kind::NODE) {
    xdmf_->writeField<IntVector_type::scalar_type>(attrs.name(), location);
    std::stringstream path;
    path << "/" << attrs.name() << "/" << cycle_;
    h5_data_->writeVector<IntVector_type::scalar_type>(path.str(), vec);
  }
}


void
OutputXDMF::write(const Teuchos::ParameterList& attrs,
                  const CompositeVector_<double_type>& vec) const
{
  // write only cell or only node entities
  if (vec.hasComponent("cell") && vec.getLocation("cell") == AmanziMesh::Entity_kind::CELL) {
    Teuchos::ParameterList cell_attrs(attrs);
    cell_attrs.set<AmanziMesh::Entity_kind>("location", AmanziMesh::Entity_kind::CELL);
    write(cell_attrs, *vec.getComponent("cell", false));
  } else if (vec.hasComponent("node") && vec.getLocation("node") == AmanziMesh::Entity_kind::NODE) {
    Teuchos::ParameterList node_attrs(attrs);
    node_attrs.set("location", AmanziMesh::Entity_kind::NODE);
    write(node_attrs, *vec.getComponent("node", false));
  }
}

void
OutputXDMF::write(const Teuchos::ParameterList& attrs, const CompositeVector_<int_type>& vec) const
{
  // write only cell or only node entities
  if (vec.hasComponent("cell") && vec.getLocation("cell") == AmanziMesh::Entity_kind::CELL) {
    Teuchos::ParameterList cell_attrs(attrs);
    cell_attrs.set("location", AmanziMesh::Entity_kind::CELL);
    write(cell_attrs, *vec.getComponent("cell", false));
  } else if (vec.hasComponent("node") && vec.getLocation("node") == AmanziMesh::Entity_kind::NODE) {
    Teuchos::ParameterList node_attrs(attrs);
    node_attrs.set("location", AmanziMesh::Entity_kind::NODE);
    write(node_attrs, *vec.getComponent("node", false));
  }
}


// can we template this (not yet...)
void
OutputXDMF::write(const Teuchos::ParameterList& attrs, const double& val) const
{
  h5_data_->writeAttribute(std::to_string(cycle_), attrs.name(), val);
}

void
OutputXDMF::write(const Teuchos::ParameterList& attrs, const int& val) const
{
  h5_data_->writeAttribute(std::to_string(cycle_), attrs.name(), val);
}

void
OutputXDMF::write(const Teuchos::ParameterList& attrs, const std::string& val) const
{
  h5_data_->writeAttribute(std::to_string(cycle_), attrs.name(), val);
}


} // namespace Amanzi
