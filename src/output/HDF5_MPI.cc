/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Erin Baker
*/

/*
  Writes data structures to HDF5 files in parallel.
*/

#include <iostream>
#include <string>

#include "HDF5_MPI.hh"

//TODO(barker): clean up debugging output
//TODO(barker): check that close file is always getting called
//TODO(barker): add error handling where appropriate
//TODO(barker): clean up formating

namespace Amanzi {

HDF5_MPI::HDF5_MPI(const Comm_ptr_type &comm, bool include_io_set)
    : viz_comm_(Teuchos::rcp_dynamic_cast<const MpiComm_type>(comm)),
      dynamic_mesh_(false),
      mesh_written_(false),
      static_mesh_cycle_(0),
      include_io_set_(include_io_set)
{
  AMANZI_ASSERT(viz_comm_.get());
  info_ = MPI_INFO_NULL;
  IOconfig_.numIOgroups = 1;
  IOconfig_.commIncoming = viz_comm_->Comm();
  parallelIO_IOgroup_init(&IOconfig_, &IOgroup_);
  mesh_ = Teuchos::null;
  NumNodes_ = 0;
  NumElems_ = 0;
  ConnLength_ = 0;
}


HDF5_MPI::HDF5_MPI(const Comm_ptr_type &comm, std::string dataFilename, bool include_io_set)
    : viz_comm_(Teuchos::rcp_dynamic_cast<const MpiComm_type>(comm)),
      H5DataFilename_(dataFilename),
      dynamic_mesh_(false),
      mesh_written_(false),
      static_mesh_cycle_(0),
      include_io_set_(include_io_set)
{
  AMANZI_ASSERT(viz_comm_.get());
  H5DataFilename_ = dataFilename;
  info_ = MPI_INFO_NULL;
  IOconfig_.numIOgroups = 1;
  IOconfig_.commIncoming = viz_comm_->Comm();
  parallelIO_IOgroup_init(&IOconfig_, &IOgroup_);
  mesh_ = Teuchos::null;
  NumNodes_ = 0;
  NumElems_ = 0;
  ConnLength_ = 0;
}


HDF5_MPI::~HDF5_MPI()
{
  parallelIO_IOgroup_cleanup(&IOgroup_);
}


void HDF5_MPI::createMeshFile(Teuchos::RCP<const AmanziMesh::Mesh> mesh, const std::string& filename)
{
  std::string xmfFilename;

  // store the mesh
  mesh_ = mesh;

  // build h5 filename
  base_filename_ = filename;
  h5Filename_ = filename + ".h5";

  // new parallel
  mesh_file_ = parallelIO_open_file(h5Filename_.c_str(), &IOgroup_, FILE_CREATE);
  if (mesh_file_ < 0) {
    Errors::Message message("HDF5_MPI::createMeshFile - error creating mesh file");
    Exceptions::amanzi_throw(message);
  }
  // close file
  parallelIO_close_file(mesh_file_, &IOgroup_);

  // Store filenames
  if (TrackXdmf() && viz_comm_->MyPID() == 0) {
    setxdmfMeshVisitFilename(filename + ".VisIt.xmf");
    // start xmf files xmlObjects stored inside functions
    createXdmfMeshVisit_();
  }
}


void HDF5_MPI::writeMesh(const double time, const int iteration)
{
  if (mesh_->is_logical()) {
    writeDualMesh(time, iteration);
    return;
  }

  const AmanziMesh::Mesh& vis_mesh = mesh_->vis_mesh();

  std::string xmfFilename;
  int globaldims[2], localdims[2];
  int *ids;

  // if this is a static mesh simulation, we only write the mesh once
  if (!dynamic_mesh_ && mesh_written_) return;

  // open the mesh file
  mesh_file_ = parallelIO_open_file(h5Filename_.c_str(), &IOgroup_, FILE_READWRITE);

  // get num_nodes, num_cells
  const Epetra_Map& nmap = vis_mesh.node_map(false);
  int nnodes_local = nmap.NumMyElements();
  int nnodes_global = nmap.NumGlobalElements();
  const Epetra_Map& ngmap = vis_mesh.node_map(true);

  const Epetra_Map& cmap = vis_mesh.cell_map(false);
  int ncells_local = cmap.NumMyElements();

  // get space dimension
  int space_dim = vis_mesh.space_dimension();
  int topo_dim = vis_mesh.manifold_dimension();

  // Get and write node coordinate info
  // -- get coords
  double *xyz = new double[nnodes_local*3];
  globaldims[0] = nnodes_global;
  globaldims[1] = 3;
  localdims[0] = nnodes_local;
  localdims[1] = 3;

  AmanziGeometry::Point xc(space_dim);
  for (int i = 0; i < nnodes_local; i++) {
    vis_mesh.node_get_coordinates(i, &xc);
    // VisIt and ParaView require all mesh entities to be in 3D space
    xyz[i*3+0] = xc[0];
    xyz[i*3+1] = xc[1];
    if (space_dim == 3) {
      xyz[i*3+2] = xc[2];
    } else {
      xyz[i*3+2] = 0.0;
    }
  }

  // -- write out coords
  std::stringstream hdf5_path;
  hdf5_path << iteration << "/Mesh/Nodes";
  // TODO(barker): add error handling: can't create/write
  parallelIO_write_dataset(xyz, PIO_DOUBLE, 2, globaldims, localdims, mesh_file_,
                           const_cast<char*>(hdf5_path.str().c_str()), &IOgroup_,
                           NONUNIFORM_CONTIGUOUS_WRITE);

  // -- clean up
  delete [] xyz;

  // -- write out node map
  ids = new int[nmap.NumMyElements()];
  for (int i=0; i<nnodes_local; i++) {
    ids[i] = nmap.GID(i);
  }
  globaldims[1] = 1;
  localdims[1] = 1;

  hdf5_path.str("");
  hdf5_path << iteration << "/Mesh/NodeMap";
  parallelIO_write_dataset(ids, PIO_INTEGER, 2, globaldims, localdims, mesh_file_,
                           const_cast<char*>(hdf5_path.str().c_str()), &IOgroup_,
                           NONUNIFORM_CONTIGUOUS_WRITE);
  delete [] ids;

  // Get and write cell-to-node connectivity information
  // -- get connectivity
  // nodes are written to h5 out of order, need info to map id to order in output
  int nnodes(nnodes_local);
  std::vector<int> nnodesAll(viz_comm_->NumProc(),0);
  viz_comm_->GatherAll(&nnodes, &nnodesAll[0], 1);
  int start(0);
  std::vector<int> startAll(viz_comm_->NumProc(),0);
  for (int i = 0; i < viz_comm_->MyPID(); i++) {
    start += nnodesAll[i];
  }
  viz_comm_->GatherAll(&start, &startAll[0],1);

  std::vector<int> gid(nnodes_global);
  std::vector<int> pid(nnodes_global);
  std::vector<int> lid(nnodes_global);
  for (int i=0; i<nnodes_global; i++) {
    gid[i] = ngmap.GID(i);
  }
  nmap.RemoteIDList(nnodes_global, &gid[0], &pid[0], &lid[0]);

  // -- determine size of connectivity vector
  // The element connectivity vector is given in the following form (XDMF spec):
  //
  //  if elem_type != POLYGON:
  //      elem_type node1 node2 ... nodeN
  //  else if space_dim == 2:
  //      POLYGON num_nodes node1 node2 ... nodeN
  //  else:
  //      // 3D POLYHEDRON, not true POLYGON
  //      // represent each element as num_faces POLYGONS, each which map to the same element ID
  //      POLYGON num_nodes_in_face1 node1 node2 ... nodeM
  //         POLYGON num_nodes_in_face2 node1 node2 ... nodeM
  //           ...
  //              POLYGON num_nodes_in_faceK node1 node2 ... nodeM
  //
  // As XDMF does not yet support POLYHEDRON, this is an attempted workaround.
  //
  // TODO(barker): make a list of cell types found,
  //               if all the same then write out as a uniform mesh of that type

  // -- pass 1: count total connections, total entities
  int local_conn(0); // length of MixedElements
  int local_entities(0); // length of ElementMap (num_cells if non-POLYHEDRON mesh)
  AmanziMesh::Entity_ID_List faces, nodes;

  for (int c=0; c!=ncells_local; ++c) {
    AmanziMesh::Cell_type ctype = vis_mesh.cell_get_type(c);
    if (getCellTypeID_(ctype) == getCellTypeID_(AmanziMesh::POLYHED)) {
      vis_mesh.cell_get_faces(c,&faces);
      for (int i=0; i!=faces.size(); ++i) {
        vis_mesh.face_get_nodes(faces[i], &nodes);
        local_conn += nodes.size() + 1;
      }
      local_conn += 2;
      local_entities++;

    } else if (getCellTypeID_(ctype) != getCellTypeID_(AmanziMesh::POLYGON)) {
      vis_mesh.cell_get_nodes(c,&nodes);
      local_conn += nodes.size() + 1;
      local_entities++;

    } else if (topo_dim == 2) {
      vis_mesh.cell_get_nodes(c,&nodes);
      local_conn += nodes.size() + 2;
      local_entities++;

    } else {
      vis_mesh.cell_get_faces(c,&faces);
      for (int i=0; i!=faces.size(); ++i) {
        vis_mesh.face_get_nodes(faces[i], &nodes);
        local_conn += nodes.size() + 2;
      }
      local_entities += faces.size();
    }
  }

  int *local_sizes = new int[2]; // length of MixedElements,
		      // length of ElementMap (num_cells if non-POLYHEDRON mesh)
  local_sizes[0] = local_conn; local_sizes[1] = local_entities;
  int *global_sizes = new int[2];
  global_sizes[0] = 0; global_sizes[1] = 0;
  viz_comm_->SumAll(local_sizes, global_sizes, 2);
  int global_conn = global_sizes[0];
  int global_entities = global_sizes[1];
  delete [] local_sizes;
  delete [] global_sizes;

  // allocate space, pass 2 to populate
  // nodeIDs need to be mapped to output IDs
  int *conn = new int[local_conn];
  int *entities = new int[local_entities];

  int lcv = 0;
  int lcv_entity = 0;
  for (int c=0; c!=ncells_local; ++c) {
    AmanziMesh::Cell_type ctype = vis_mesh.cell_get_type(c);
    if (getCellTypeID_(ctype) == getCellTypeID_(AmanziMesh::POLYHED)) {
      vis_mesh.cell_get_faces(c,&faces);

      // store cell type id and number of faces
      conn[lcv++] = getCellTypeID_(ctype);
      conn[lcv++] = faces.size();

      for (int j=0; j!=faces.size(); ++j) {
        // store node count, then nodes in the correct order
        vis_mesh.face_get_nodes(faces[j], &nodes);
        conn[lcv++] = nodes.size();

        for (int i=0; i!=nodes.size(); ++i) {
          if (nmap.MyLID(nodes[i])) {
            conn[lcv++] = nodes[i] + startAll[viz_comm_->MyPID()];
          } else {
            conn[lcv++] = lid[nodes[i]] + startAll[pid[nodes[i]]];
          }
        }
      }

      // store entity
      entities[lcv_entity++] = cmap.GID(c);

    } else if (getCellTypeID_(ctype) != getCellTypeID_(AmanziMesh::POLYGON)) {
      // store cell type id
      conn[lcv++] = getCellTypeID_(ctype);

      // store nodes in the correct order
      vis_mesh.cell_get_nodes(c, &nodes);

      for (int i=0; i!=nodes.size(); ++i) {
	if (nmap.MyLID(nodes[i])) {
	  conn[lcv++] = nodes[i] + startAll[viz_comm_->MyPID()];
	} else {
	  conn[lcv++] = lid[nodes[i]] + startAll[pid[nodes[i]]];
	}
      }

      // store entity
      entities[lcv_entity++] = cmap.GID(c);
      
    } else if (topo_dim == 2) {
      // store cell type id
      conn[lcv++] = getCellTypeID_(ctype);

      // store node count, then nodes in the correct order
      vis_mesh.cell_get_nodes(c, &nodes);
      conn[lcv++] = nodes.size();

      for (int i=0; i!=nodes.size(); ++i) {
	if (nmap.MyLID(nodes[i])) {
	  conn[lcv++] = nodes[i] + startAll[viz_comm_->MyPID()];
	} else {
	  conn[lcv++] = lid[nodes[i]] + startAll[pid[nodes[i]]];
	}
      }

      // store entity
      entities[lcv_entity++] = cmap.GID(c);
      
    } else {
      vis_mesh.cell_get_faces(c,&faces);

      for (int j=0; j!=faces.size(); ++j) {
	// store cell type id
	conn[lcv++] = getCellTypeID_(ctype);

	// store node count, then nodes in the correct order
	vis_mesh.face_get_nodes(faces[j], &nodes);
	conn[lcv++] = nodes.size();

	for (int i=0; i!=nodes.size(); ++i) {
	  if (nmap.MyLID(nodes[i])) {
	    conn[lcv++] = nodes[i] + startAll[viz_comm_->MyPID()];
	  } else {
	    conn[lcv++] = lid[nodes[i]] + startAll[pid[nodes[i]]];
	  }
	}

	// store entity
	entities[lcv_entity++] = cmap.GID(c);
      }
    }
  }
  
  // write out connectivity
  hdf5_path.str("");
  hdf5_path << iteration << "/Mesh/MixedElements";
  globaldims[0] = global_conn; globaldims[1] = 1;
  localdims[0] = local_conn; localdims[1] = 1;
  parallelIO_write_dataset(conn, PIO_INTEGER, 2, globaldims, localdims, mesh_file_,
                           const_cast<char*>(hdf5_path.str().c_str()), &IOgroup_,
                           NONUNIFORM_CONTIGUOUS_WRITE);
  delete [] conn;
  hdf5_path.flush();

  // write out entity map
  hdf5_path.str("");
  hdf5_path << iteration << "/Mesh/ElementMap";
  globaldims[0] = global_entities; globaldims[1] = 1;
  localdims[0] = local_entities; localdims[1] = 1;
  parallelIO_write_dataset(entities, PIO_INTEGER, 2, globaldims, localdims, mesh_file_,
                           const_cast<char*>(hdf5_path.str().c_str()), &IOgroup_,
                           NONUNIFORM_CONTIGUOUS_WRITE);
  delete [] entities;

  // close file
  parallelIO_close_file(mesh_file_, &IOgroup_);

  // Store information
  setH5MeshFilename(h5Filename_);
  setNumNodes(nnodes_global);
  setNumElems(global_entities);
  setConnLength(global_conn);

  // Create and write out accompanying Xdmf file
  if (TrackXdmf() && viz_comm_->MyPID() == 0) {
    //TODO(barker): if implement type tracking, then update this as needed

    // This can be optimized in the case where all elements are of the same type
    cname_ = "Mixed";
    xmfFilename = base_filename_ + ".xmf";
    createXdmfMesh_(base_filename_, time, iteration);

    // update Mesh VisIt xdmf files
    std::stringstream fname;
    fname << base_filename_ << ".h5." << iteration << ".xmf";

    writeXdmfMeshVisitGrid_(fname.str());
    std::ofstream of;
    of.open(xdmfMeshVisitFilename().c_str());
    of << xmlMeshVisit();
    of.close();
  }

  mesh_written_ = true;
  static_mesh_cycle_ = iteration;
}


void HDF5_MPI::writeDualMesh(const double time, const int iteration)
{
  std::string xmfFilename;
  int globaldims[2], localdims[2];
  int *ids;

  const AmanziMesh::Mesh& vis_mesh = mesh_->vis_mesh();

  // if this is a static mesh simulation, we only write the mesh once
  if (!dynamic_mesh_ && mesh_written_) return;

  // open the mesh file
  mesh_file_ = parallelIO_open_file(h5Filename_.c_str(), &IOgroup_, FILE_READWRITE);

  // get num_nodes, num_cells
  const Epetra_Map &nmap = vis_mesh.cell_map(false);
  int nnodes_local = nmap.NumMyElements();
  int nnodes_global = nmap.NumGlobalElements();
  const Epetra_Map &ngmap = vis_mesh.cell_map(true);

  // get space dimension
  int space_dim = vis_mesh.space_dimension();

  // Get and write node coordinate info
  // -- get coords
  double *nodes = new double[nnodes_local*3];
  globaldims[0] = nnodes_global;
  globaldims[1] = 3;
  localdims[0] = nnodes_local;
  localdims[1] = 3;

  for (int i = 0; i < nnodes_local; i++) {
    const AmanziGeometry::Point& xc = vis_mesh.cell_centroid(i);
    // VisIt and ParaView require all mesh entities to be in 3D space
    nodes[i*3+0] = xc[0];
    nodes[i*3+1] = xc[1];
    if (space_dim == 3) {
      nodes[i*3+2] = xc[2];
    } else {
      nodes[i*3+2] = 0.0;
    }
  }

  // -- write out coords
  std::stringstream hdf5_path;
  hdf5_path << iteration << "/Mesh/Nodes";
  // TODO(barker): add error handling: can't create/write
  parallelIO_write_dataset(nodes, PIO_DOUBLE, 2, globaldims, localdims, mesh_file_,
                           const_cast<char*>(hdf5_path.str().c_str()), &IOgroup_,
                           NONUNIFORM_CONTIGUOUS_WRITE);

  // -- clean up
  delete [] nodes;

  // -- write out node map
  ids = new int[nmap.NumMyElements()];
  for (int i=0; i<nnodes_local; i++) {
    ids[i] = nmap.GID(i);
  }
  globaldims[1] = 1;
  localdims[1] = 1;

  hdf5_path.str("");
  hdf5_path << iteration << "/Mesh/NodeMap";
  parallelIO_write_dataset(ids, PIO_INTEGER, 2, globaldims, localdims, mesh_file_,
                           const_cast<char*>(hdf5_path.str().c_str()), &IOgroup_,
                           NONUNIFORM_CONTIGUOUS_WRITE);
  delete [] ids;

  // Get and write cell-to-node connectivity information
  // -- get connectivity
  // nodes are written to h5 out of order, need info to map id to order in output
  int nnodes(nnodes_local);
  std::vector<int> nnodesAll(viz_comm_->NumProc(),0);
  viz_comm_->GatherAll(&nnodes, &nnodesAll[0], 1);
  int start(0);
  std::vector<int> startAll(viz_comm_->NumProc(),0);
  for (int i = 0; i < viz_comm_->MyPID(); i++) {
    start += nnodesAll[i];
  }
  viz_comm_->GatherAll(&start, &startAll[0],1);

  std::vector<int> gid(nnodes_global);
  std::vector<int> pid(nnodes_global);
  std::vector<int> lid(nnodes_global);
  for (int i=0; i<nnodes_global; i++) {
    gid[i] = ngmap.GID(i);
  }
  nmap.RemoteIDList(nnodes_global, &gid[0], &pid[0], &lid[0]);

  // -- pass 1: count total connections, total entities
  // For the dual, connections are all faces with > 1 cells.
  int local_conn(0); // length of MixedElements
  int local_entities(0); // length of ElementMap (num_cells if non-POLYHEDRON mesh)
  AmanziMesh::Entity_ID_List cells;

  const Epetra_Map &fmap = vis_mesh.face_map(false);
  int nfaces_local = fmap.NumMyElements();

  for (int f=0; f!=nfaces_local; ++f) {
    vis_mesh.face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
    if (cells.size() > 1) {
      local_conn += cells.size();
      local_entities++;
    }
  }

  int *local_sizes = new int[2]; // length of MixedElements,
		      // length of ElementMap (num_cells if non-POLYHEDRON mesh)
  local_sizes[0] = local_conn; local_sizes[1] = local_entities;
  int *global_sizes = new int[2];
  global_sizes[0] = 0; global_sizes[1] = 0;
  viz_comm_->SumAll(local_sizes, global_sizes, 2);
  int global_conn = global_sizes[0];
  int global_entities = global_sizes[1];
  delete [] local_sizes;
  delete [] global_sizes;

  // allocate space, pass 2 to populate
  // nodeIDs need to be mapped to output IDs
  int *conn = new int[local_conn];
  int *entities = new int[local_entities];

  int lcv = 0;
  int lcv_entity = 0;
  int internal_f = 0;
  for (int f=0; f!=nfaces_local; ++f) {
    vis_mesh.face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);
    if (cells.size() > 1) {
      // store cell type id
      // conn[lcv++] = 2;

      // store nodes in the correct order
      for (int i=0; i!=cells.size(); ++i) {
        if (nmap.MyLID(cells[i])) {
          conn[lcv++] = cells[i] + startAll[viz_comm_->MyPID()];
        } else {
          conn[lcv++] = lid[cells[i]] + startAll[pid[cells[i]]];
        }
      }

      // store entity
      entities[lcv_entity++] = startAll[viz_comm_->MyPID()] + internal_f;
      internal_f++;
      
    // } else if (space_dim == 2) {
    //   // store cell type id
    //   conn[lcv++] = getCellTypeID_(ctype);

    //   // store node count, then nodes in the correct order
    //   AmanziMesh::Entity_ID_List nodes;
    //   vis_mesh.cell_get_nodes(c, &nodes);
    //   conn[lcv++] = nodes.size();

    //   for (int i=0; i!=nodes.size(); ++i) {
    //     if (nmap.MyLID(nodes[i])) {
    //       conn[lcv++] = nodes[i] + startAll[viz_comm_->MyPID()];
    //     } else {
    //       conn[lcv++] = lid[nodes[i]] + startAll[pid[nodes[i]]];
    //     }
    //   }

    //   // store entity
    //   entities[lcv_entity++] = cmap.GID(c);
    }
  }
  
  // write out connectivity
  hdf5_path.str("");
  hdf5_path << iteration << "/Mesh/MixedElements";
  globaldims[0] = global_conn; globaldims[1] = 1;
  localdims[0] = local_conn; localdims[1] = 1;
  parallelIO_write_dataset(conn, PIO_INTEGER, 2, globaldims, localdims, mesh_file_,
                           const_cast<char*>(hdf5_path.str().c_str()), &IOgroup_,
                           NONUNIFORM_CONTIGUOUS_WRITE);
  delete [] conn;
  hdf5_path.flush();

  // write out entity map
  hdf5_path.str("");
  hdf5_path << iteration << "/Mesh/ElementMap";
  globaldims[0] = global_entities; globaldims[1] = 1;
  localdims[0] = local_entities; localdims[1] = 1;
  parallelIO_write_dataset(entities, PIO_INTEGER, 2, globaldims, localdims, mesh_file_,
                           const_cast<char*>(hdf5_path.str().c_str()), &IOgroup_,
                           NONUNIFORM_CONTIGUOUS_WRITE);
  delete [] entities;

  // close file
  parallelIO_close_file(mesh_file_, &IOgroup_);

  // Store information
  setH5MeshFilename(h5Filename_);
  setNumNodes(nnodes_global);
  setNumElems(global_entities);
  setConnLength(global_conn);

  // Create and write out accompanying Xdmf file
  if (TrackXdmf() && viz_comm_->MyPID() == 0) {
    //TODO(barker): if implement type tracking, then update this as needed

    // This can be optimized in the case where all elements are of the same type
    cname_ = "Polyline";
    xmfFilename = base_filename_ + ".xmf";
    createXdmfMesh_(base_filename_, time, iteration);

    // update Mesh VisIt xdmf files
    std::stringstream fname;
    fname << base_filename_ << ".h5." << iteration << ".xmf";

    writeXdmfMeshVisitGrid_(fname.str());
    std::ofstream of;
    of.open(xdmfMeshVisitFilename().c_str());
    of << xmlMeshVisit();
    of.close();
  }

  mesh_written_ = true;
  static_mesh_cycle_ = iteration;
}


void HDF5_MPI::createDataFile(const std::string& base_filename)
{
  // ?? input mesh filename or grab global mesh filename
  // ->assumes global name exists!!
  // build h5 filename
  std::string h5filename = base_filename + ".h5";

  // new parallel
  data_file_ = parallelIO_open_file(h5filename.c_str(), &IOgroup_, FILE_CREATE);

  if (data_file_ < 0) {
    Errors::Message message("HDF5_MPI::createDataFile - error creating data file");
    Exceptions::amanzi_throw(message);
  }

  // close file
  parallelIO_close_file(data_file_, &IOgroup_);

  // Store filenames
  setH5DataFilename(h5filename);
  if (TrackXdmf() && viz_comm_->MyPID() == 0) {
    xdmfVisitFilename_ = base_filename + ".VisIt.xmf";
    xdmfVisitBaseFilename_ = base_filename;
    // start xmf files xmlObjects stored inside functions
    createXdmfVisit_();
  }
}


void HDF5_MPI::open_h5file(bool read_only)
{
  if (read_only) {
    data_file_ = parallelIO_open_file(H5DataFilename_.c_str(), &IOgroup_,
             FILE_READWRITE);
    if (data_file_ < 0) {
      Errors::Message msg;
      msg << "HDF5_MPI: error opening file \"" << H5DataFilename_ << "\" with READ_WRITE access.";
      Exceptions::amanzi_throw(msg);
    }
  } else {
    data_file_ = parallelIO_open_file(H5DataFilename_.c_str(), &IOgroup_,
             FILE_READWRITE);
    if (data_file_ < 0) {
      Errors::Message msg;
      msg << "HDF5_MPI: error opening file \"" << H5DataFilename_ << "\" with READ_ONLY access.";
      Exceptions::amanzi_throw(msg);
    }
  }
}


void HDF5_MPI::close_h5file() {
  parallelIO_close_file(data_file_, &IOgroup_);
  data_file_ = -1;
}


void HDF5_MPI::createTimestep(double time, int iteration, const std::string& tag)
{
  setIteration(iteration);
  setTime(time);
  set_tag(tag);

  if (TrackXdmf() && viz_comm_->MyPID() == 0) {
    // create single step xdmf file
    Teuchos::XMLObject tmp("Xdmf");
    tmp.addChild(addXdmfHeaderLocal_("Mesh",time,iteration));
    std::stringstream filename;
    filename << H5DataFilename() << "." << iteration << ".xmf";
    of_timestep_.open(filename.str().c_str());
    // channel will be closed when the endTimestep() is called
    setxdmfStepFilename(filename.str());

    // Store information
    xmlStep_ = tmp;
  }
}


void HDF5_MPI::endTimestep()
{
  if (TrackXdmf() && viz_comm_->MyPID() == 0) {
    of_timestep_ << xmlStep_;
    of_timestep_.close();

    // compare ifields in two subsequent steps and open a new meta-data 
    // channel if they differ
    auto fields = extractFields_(xmlStep_);
    auto fields_prev = extractFields_(xmlStep_prev_);
    if (fields != fields_prev && !xmlStep_prev_.isEmpty()) {
      createXdmfVisit_();
    }

    xmlStep_prev_ = xmlStep_;

    // add a new time step to global VisIt xdmf files
    // TODO(barker): how to get to grid collection node, rather than root???
    std::string record = H5DataFilename() + "." + std::to_string(Iteration()) + ".xmf";
    writeXdmfVisitGrid_(record);
    // TODO(barker): where to write out depends on where the root node is
    // ?? how to terminate stream or switch to new file out??
    int count(0);
    for (auto& xmf : xmlVisit_) {
      std::string full_name;
      if (include_io_set_) {
        full_name = xdmfVisitBaseFilename_ + "_io-set" + std::to_string(count++) + ".VisIt.xmf";
      } else {
        full_name = xdmfVisitBaseFilename_ + ".VisIt.xmf";
      }
      std::ofstream of;
      of.open(full_name.c_str());
      of << xmf;
      of.close();
    }
  }
}


std::set<std::string> HDF5_MPI::extractFields_(const Teuchos::XMLObject& xml)
{
  std::set<std::string> fields;

  if (!xml.isEmpty()) {
    auto node = findMeshNode_(xml);
    for (int i = 0; i < node.numChildren(); i++) {
      auto tmp = node.getChild(i);
      if (tmp.getTag() == "Attribute" && tmp.hasAttribute("Name") && tmp.hasAttribute("Type")) {
        if (tmp.getAttribute("Type") == "Scalar") {
          fields.insert(tmp.getAttribute("Name"));
        }
      }
    }
  }
  return fields;
}


void HDF5_MPI::writeAttrString(const std::string value, const std::string attrname)
{
  std::string h5path = "/";

  char *loc_attrname = new char[attrname.size()+1];
  strcpy(loc_attrname,attrname.c_str());

  char *loc_h5path = new char[h5path.size()+1];
  strcpy(loc_h5path,h5path.c_str());

  char *loc_value = new char[value.size()+1];
  strcpy(loc_value, value.c_str());


  parallelIO_write_simple_attr(loc_attrname,
                               loc_value,
                               PIO_STRING,
                               data_file_,
                               loc_h5path,
                               &IOgroup_);
  delete [] loc_value;
  delete [] loc_h5path;
  delete [] loc_attrname;
}


void HDF5_MPI::writeAttrReal(double value, const std::string attrname)
{
  std::string h5path = "/";

  char *loc_attrname = new char[attrname.size()+1];
  strcpy(loc_attrname,attrname.c_str());

  char *loc_h5path = new char[h5path.size()+1];
  strcpy(loc_h5path,h5path.c_str());

  parallelIO_write_simple_attr(loc_attrname,
                               &value,
                               PIO_DOUBLE,
                               data_file_,
                               loc_h5path,
                               &IOgroup_);
  delete [] loc_h5path;
  delete [] loc_attrname;
}


void HDF5_MPI::writeAttrReal(double value, const std::string attrname, std::string h5path)
{
  char *loc_attrname = new char[attrname.size()+1];
  strcpy(loc_attrname,attrname.c_str());

  char *loc_h5path = new char[h5path.size()+1];
  strcpy(loc_h5path,h5path.c_str());

  parallelIO_write_simple_attr(loc_attrname,
                               &value,
                               PIO_DOUBLE,
                               data_file_,
                               loc_h5path,
                               &IOgroup_);
  delete [] loc_h5path;
  delete [] loc_attrname;
}


void HDF5_MPI::writeAttrReal(double* value, int ndim, const std::string attrname)
{
  std::string h5path = "/";

  char *loc_attrname = new char[attrname.size()+1];
  strcpy(loc_attrname,attrname.c_str());

  char *loc_h5path = new char[h5path.size()+1];
  strcpy(loc_h5path,h5path.c_str());

  int ndims(1);
  int *adims = new int[1];
  adims[0] = ndim;

  parallelIO_write_attr(loc_attrname,
                        reinterpret_cast<void*>(value),
                        PIO_DOUBLE,
                        ndims,
                        adims,
                        data_file_,
                        loc_h5path,
                        &IOgroup_);

  delete [] loc_h5path;
  delete [] loc_attrname;
  delete [] adims;
}


void HDF5_MPI::writeAttrInt(int value, const std::string attrname)
{
  std::string h5path = "/";

  char *loc_attrname = new char[attrname.size()+1];
  strcpy(loc_attrname,attrname.c_str());

  char *loc_h5path = new char[h5path.size()+1];
  strcpy(loc_h5path,h5path.c_str());

  parallelIO_write_simple_attr(loc_attrname,
                               &value,
                               PIO_INTEGER,
                               data_file_,
                               loc_h5path,
                               &IOgroup_);
  delete [] loc_h5path;
  delete [] loc_attrname;
}


void HDF5_MPI::writeAttrInt(int* value, int ndim, const std::string attrname)
{
  std::string h5path = "/";

  char *loc_attrname = new char[attrname.size()+1];
  strcpy(loc_attrname,attrname.c_str());

  char *loc_h5path = new char[h5path.size()+1];
  strcpy(loc_h5path,h5path.c_str());

  int ndims(1);
  int *adims = new int[1];
  adims[0] = ndim;

  parallelIO_write_attr(loc_attrname,
                        reinterpret_cast<void*>(value),
                        PIO_INTEGER,
                        ndims,
                        adims,
                        data_file_,
                        loc_h5path,
                        &IOgroup_);

  delete [] loc_h5path;
  delete [] loc_attrname;
  delete [] adims;
}


void HDF5_MPI::readAttrString(std::string &value, const std::string attrname)
{
  std::string h5path = "/";

  char *loc_attrname = new char[attrname.size()+1];
  strcpy(loc_attrname,attrname.c_str());

  char *loc_h5path = new char[h5path.size()+1];
  strcpy(loc_h5path,h5path.c_str());

  char *loc_value;
  
  parallelIO_read_simple_attr(loc_attrname,
                              reinterpret_cast<void**>(&loc_value),
                              PIO_STRING,
                              data_file_,
                              loc_h5path,
                              &IOgroup_);

  value = std::string(loc_value);

  free(loc_value);
  delete [] loc_h5path;
  delete [] loc_attrname;
}


void HDF5_MPI::readAttrReal(double &value, const std::string attrname)
{
  std::string h5path = "/";

  char *loc_attrname = new char[attrname.size()+1];
  strcpy(loc_attrname,attrname.c_str());

  char *loc_h5path = new char[h5path.size()+1];
  strcpy(loc_h5path,h5path.c_str());

  double *loc_value;
  
  parallelIO_read_simple_attr(loc_attrname,
                              reinterpret_cast<void**>(&loc_value),
                              PIO_DOUBLE,
                              data_file_,
                              loc_h5path,
                              &IOgroup_);

  value = *loc_value;

  free(loc_value);
  delete [] loc_h5path;
  delete [] loc_attrname;
}
  

void HDF5_MPI::readAttrReal(double **value, int *ndim, const std::string attrname)
{
  std::string h5path = "/";

  char *loc_attrname = new char[attrname.size()+1];
  strcpy(loc_attrname,attrname.c_str());

  char *loc_h5path = new char[h5path.size()+1];
  strcpy(loc_h5path,h5path.c_str());

  double *loc_value;
  int *pdims;
  int ndims;
  
  parallelIO_read_attr(loc_attrname,
                       reinterpret_cast<void**>(&loc_value),
                       PIO_DOUBLE,
                       &ndims,
                       &pdims,
                       data_file_,
                       loc_h5path,
                       &IOgroup_);

  *value = loc_value;
  *ndim = pdims[0];  // works only for one-dimensional vectors.

  free(pdims);
  delete [] loc_h5path;
  delete [] loc_attrname;
}


void HDF5_MPI::readAttrInt(int &value, const std::string attrname)
{
  std::string h5path = "/";

  char *loc_attrname = new char[attrname.size()+1];
  strcpy(loc_attrname,attrname.c_str());

  char *loc_h5path = new char[h5path.size()+1];
  strcpy(loc_h5path,h5path.c_str());

  int *loc_value;
  
  parallelIO_read_simple_attr(loc_attrname,
                              reinterpret_cast<void**>(&loc_value),
                              PIO_INTEGER,
                              data_file_,
                              loc_h5path,
                              &IOgroup_);

  value = *loc_value;

  free(loc_value);
  delete [] loc_h5path;
  delete [] loc_attrname;
}


void HDF5_MPI::readAttrInt(int **value, int *ndim, const std::string attrname)
{
  std::string h5path = "/";

  char *loc_attrname = new char[attrname.size()+1];
  strcpy(loc_attrname,attrname.c_str());

  char *loc_h5path = new char[h5path.size()+1];
  strcpy(loc_h5path,h5path.c_str());

  int *loc_value, *pdims;
  int ndims;
  
  parallelIO_read_attr(loc_attrname,
                       reinterpret_cast<void**>(&loc_value),
                       PIO_INTEGER,
                       &ndims,
                       &pdims,
                       data_file_,
                       loc_h5path,
                       &IOgroup_);

  *value = loc_value;
  *ndim = pdims[0];  // works only for one-dimensional vectors.

  free(pdims);
  delete [] loc_h5path;
  delete [] loc_attrname;
}


void HDF5_MPI::writeDataString(char **x, int num_entries, const std::string& varname)
{
  char *h5path = new char[varname.size() + 1];
  strcpy(h5path, varname.c_str());

  parallelIO_write_str_array(x, num_entries, data_file_, h5path, &IOgroup_);

  delete [] h5path;
}


void HDF5_MPI::readDataString(char ***x, int *num_entries, const std::string& varname)
{
  char *h5path = new char[varname.size() + 1];
  strcpy(h5path, varname.c_str());
  int ndims, dims[2], tmpsize;

  bool exists = checkFieldData_(h5path);
  if (exists) {
    parallelIO_get_dataset_ndims(&ndims, data_file_, h5path, &IOgroup_);
    parallelIO_get_dataset_dims(dims, data_file_, h5path, &IOgroup_);
    parallelIO_get_dataset_size(&tmpsize, data_file_, h5path, &IOgroup_);

    char **strData;
    parallelIO_read_str_array(&strData, &tmpsize, data_file_, h5path, &IOgroup_);

    *x = strData;
    *num_entries = tmpsize;
  } else {
    *num_entries = 0;
  }

  delete [] h5path;
}


void HDF5_MPI::writeDataReal(const Epetra_Vector &x, const std::string& varname) {
  writeFieldData_(x, varname, PIO_DOUBLE, "NONE");
}


void HDF5_MPI::writeDataInt(const Epetra_Vector &x, const std::string& varname) {
  writeFieldData_(x, varname, PIO_INTEGER, "NONE");
}


void HDF5_MPI::writeCellDataReal(const Epetra_Vector &x, const std::string& varname) {
  writeFieldData_(x, varname, PIO_DOUBLE, "Cell");
}


void HDF5_MPI::writeCellDataInt(const Epetra_Vector &x, const std::string& varname) {
  writeFieldData_(x, varname, PIO_INTEGER, "Cell");
}


void HDF5_MPI::writeNodeDataReal(const Epetra_Vector &x, const std::string& varname) {
  writeFieldData_(x, varname, PIO_DOUBLE, "Node");
}


void HDF5_MPI::writeNodeDataInt(const Epetra_Vector &x, const std::string& varname) {
  writeFieldData_(x, varname, PIO_INTEGER, "Node");
}


void HDF5_MPI::writeFieldData_(const Epetra_Vector &x, const std::string& varname,
                               datatype_t type, std::string loc)
{
  // write field data
  double *data;
  x.ExtractView(&data);

  int globaldims[2];
  int localdims[2];
  globaldims[0] = x.GlobalLength();
  globaldims[1] = 1;
  localdims[0] = x.MyLength();
  localdims[1] = 1;

  // TODO(barker): how to build path name?? probably still need iteration number
  std::stringstream h5path;
  h5path << varname;

  // TODO(barker): add error handling: can't write/create

  //hid_t file = parallelIO_open_file(H5DataFilename_.c_str(), &IOgroup_, FILE_READWRITE);
  //MB: if (file < 0) {
  //MB:   Errors::Message message("HDF5_MPI::writeFieldData_ - error opening data file to write field data");
  //MB:   Exceptions::amanzi_throw(message);
  //MB: }

  if (TrackXdmf()) {
    h5path << "/" << Iteration();
  }

  char *tmp;
  tmp = new char[h5path.str().size() + 1];
  strcpy(tmp,h5path.str().c_str());

  parallelIO_write_dataset(data, type, 2, globaldims, localdims, data_file_, tmp,
                           &IOgroup_, NONUNIFORM_CONTIGUOUS_WRITE);

  // TODO(barker): add error handling: can't write
  if (TrackXdmf() ) {
    writeAttrReal(Time(), "Time", h5path.str());
    if (viz_comm_->MyPID() == 0) {
      // TODO(barker): get grid node, node.addChild(addXdmfAttribute)
      Teuchos::XMLObject node = findMeshNode_(xmlStep_);
      node.addChild(addXdmfAttribute_(varname, loc, globaldims[0], h5path.str()));
    }
  }
  delete [] tmp;
}


void HDF5_MPI::writeDatasetReal(double* data, int nloc, int nglb, const std::string& varname)
{
  int globaldims[2], localdims[2];
  globaldims[0] = nglb;
  globaldims[1] = 1;
  localdims[0] = nloc;
  localdims[1] = 1;

  char *tmp;
  tmp = new char[varname.size() + 1];
  strcpy(tmp, varname.c_str());

  parallelIO_write_dataset(data, PIO_DOUBLE, 2, globaldims, localdims, data_file_, tmp,
                           &IOgroup_, NONUNIFORM_CONTIGUOUS_WRITE);

  delete [] tmp;
}


bool HDF5_MPI::readData(Epetra_Vector &x, const std::string varname) {
  return readFieldData_(x, varname, PIO_DOUBLE);
}


bool HDF5_MPI::checkFieldData_(const std::string& varname)
{
  char *h5path = new char[varname.size() + 1];
  strcpy(h5path, varname.c_str());
  bool exists = false;

  if (viz_comm_->MyPID() != 0) {
    MPI_Bcast(&exists, 1, MPI_C_BOOL, 0, viz_comm_->Comm());
  } else {
    iofile_t *currfile;
    currfile = IOgroup_.file[data_file_];
    // old interface has issues with CLang 10.0.1
    // exists = H5Lexists(currfile->fid, h5path, H5P_DEFAULT);
    exists = parallelIO_name_exists(currfile->fid, h5path);

    if (!exists) {
      std::cout << "Field " << h5path << " is not found in hdf5 file.\n";
    }

    MPI_Bcast(&exists, 1, MPI_C_BOOL, 0, viz_comm_->Comm()); 
  } 
  
  delete[] h5path;

  return exists;
}


bool HDF5_MPI::readFieldData_(Epetra_Vector &x, const std::string& varname,
                              datatype_t type)
{
  if (!checkFieldData_(varname)) return false;

  char *h5path = new char[varname.size() + 1];
  strcpy(h5path, varname.c_str());

  int ndims;
  parallelIO_get_dataset_ndims(&ndims, data_file_, h5path, &IOgroup_);
  
  if (ndims < 0) {
    if (viz_comm_->MyPID() == 0) {
      std::cout<< "Dimension of the field "<<h5path<<" is negative.\n";
    }
    return false;
  }

  int globaldims[ndims], localdims[ndims];
  parallelIO_get_dataset_dims(globaldims, data_file_, h5path, &IOgroup_);
  localdims[0] = x.MyLength();
  localdims[1] = globaldims[1];

  double *data = new double[localdims[0]*localdims[1]];
  parallelIO_read_dataset(data, type, ndims, globaldims, localdims,
                          data_file_, h5path, &IOgroup_, NONUNIFORM_CONTIGUOUS_READ);

  // Trilinos' ReplaceMyValues() works with elements only and cannot 
  // be used here for points
  for (int i=0; i<localdims[0]; ++i) x[i] = data[i];

  delete [] data;
  delete [] h5path;

  return true;
}


bool HDF5_MPI::readDatasetReal(double **data, int nloc, const std::string& varname)
{
  char *h5path = new char[varname.size() + 1];
  strcpy(h5path, varname.c_str());

  if (!checkFieldData_(varname)) return false;

  int ndims;
  parallelIO_get_dataset_ndims(&ndims, data_file_, h5path, &IOgroup_);
  
  if (ndims < 0) {
    if (viz_comm_->MyPID() == 0) {
      std::cout<< "Dimension of the dataset "<<h5path<<" is negative.\n";
    }
    return false;
  }

  // root will read all data
  int globaldims[ndims], localdims[ndims];
  parallelIO_get_dataset_dims(globaldims, data_file_, h5path, &IOgroup_);
  localdims[0] = nloc;
  localdims[1] = globaldims[1];

  int size = std::max(1, localdims[0] * localdims[1]);
  *data = new double[size];
  parallelIO_read_dataset(*data, PIO_DOUBLE, ndims, globaldims, localdims,
                          data_file_, h5path, &IOgroup_, NONUNIFORM_CONTIGUOUS_READ);

  delete [] h5path;

  return true;
}


int HDF5_MPI::getCellTypeID_(AmanziMesh::Cell_type type)
{
  //TODO(barker): how to return polyhedra?
  // cell type id's defined in Xdmf/include/XdmfTopology.h

  AMANZI_ASSERT (cell_valid_type(type));

  switch (type)
  {
    case AmanziMesh::POLYGON:
      return 3;
    case AmanziMesh::TRI:
      return 4;
    case AmanziMesh::QUAD:
      return 5;
    case AmanziMesh::TET:
      return 6;
    case AmanziMesh::PYRAMID:
      return 7;
    case AmanziMesh::PRISM:
      return 8; //wedge
    case AmanziMesh::HEX:
      return 9;
    case AmanziMesh::POLYHED:
      return 16; // see http://www.xdmf.org/index.php/XDMF_Model_and_Format
    default:
      return 3; // unknown, for now same as polygon
  }
}


void HDF5_MPI::createXdmfMesh_(const std::string filename,
                               const double time, const int iteration)
{
  // TODO(barker): add error handling: can't open/write
  Teuchos::XMLObject mesh("Xdmf");

  std::stringstream mesh_name;
  mesh_name << "Mesh " << iteration;
  std::stringstream fname;
  fname << filename << ".h5." << iteration << ".xmf";

  // build xml object
  mesh.addChild(addXdmfHeaderLocal_(mesh_name.str().c_str(),time,iteration));

  // write xmf
  std::ofstream of(fname.str().c_str());
  of << HDF5_MPI::xdmfHeader_ << mesh << std::endl;
  of.close();
}


void HDF5_MPI::createXdmfMeshVisit_()
{
  Teuchos::XMLObject xmf("Xdmf");
  xmf.addAttribute("xmlns:xi", "http://www.w3.org/2001/XInclude");
  xmf.addAttribute("Version", "2.0");

  // build xml object
  xmf.addChild(addXdmfHeaderGlobal_());

  // write xmf
  std::ofstream of(xdmfMeshVisitFilename().c_str());
  of << HDF5_MPI::xdmfHeader_ << xmf << std::endl;
  of.close();

  // Store VisIt XMLObject
  xmlMeshVisit_ = xmf;
}


void HDF5_MPI::createXdmfVisit_()
{
  // TODO(barker): add error handling: can't open/write
  Teuchos::XMLObject xmf("Xdmf");
  xmf.addAttribute("xmlns:xi", "http://www.w3.org/2001/XInclude");
  xmf.addAttribute("Version", "2.0");

  // build xml object
  xmf.addChild(addXdmfHeaderGlobal_());

  // write xmf
  std::string full_name;
  if (include_io_set_) {
    int count = xmlVisit_.size();
    full_name = xdmfVisitBaseFilename_ + "_io-set" + std::to_string(count++) + ".VisIt.xmf";
  } else {
    full_name = xdmfVisitBaseFilename_ + ".VisIt.xmf";
  }
  std::ofstream of(full_name.c_str());
  of << HDF5_MPI::xdmfHeader_ << xmf << std::endl;
  of.close();

  // Store VisIt XMLObject
  xmlVisit_.push_back(xmf);
}


Teuchos::XMLObject HDF5_MPI::addXdmfHeaderGlobal_()
{
  Teuchos::XMLObject domain("Domain");

  Teuchos::XMLObject grid("Grid");
  grid.addAttribute("GridType", "Collection");
  grid.addAttribute("CollectionType", "Temporal");
  domain.addChild(grid);

  return domain;
}


Teuchos::XMLObject HDF5_MPI::addXdmfHeaderLocal_(
    const std::string& name, const double value, const int cycle)
{
  Teuchos::XMLObject domain("Domain");

  Teuchos::XMLObject grid("Grid");
  grid.addAttribute("Name", name);
  domain.addChild(grid);
  grid.addChild(addXdmfTopo_(cycle));
  grid.addChild(addXdmfGeo_(cycle));

  Teuchos::XMLObject time("Time");
  time.addDouble("Value", value);
  grid.addChild(time);

  return domain;
}


Teuchos::XMLObject HDF5_MPI::addXdmfTopo_(const int cycle)
{
  std::stringstream tmp, tmp1;

  // NEW MIXED MESH
  //TODO(barker): need to pass in topotype - or assume Mixed always
  //TODO(barker): need to pass in connectivity length, or store somewhere
  Teuchos::XMLObject topo("Topology");
  topo.addAttribute("TopologyType", cname_);
  topo.addInt("NumberOfElements", NumElems());
  topo.addAttribute("Name", "mixedtopo");

  Teuchos::XMLObject DataItem("DataItem");
  DataItem.addAttribute("DataType", "Int");
  DataItem.addInt("Dimensions", ConnLength());
  DataItem.addAttribute("Format", "HDF");

  if (dynamic_mesh_) {
    tmp1 << stripFilename_(H5MeshFilename()) << ":/" << cycle << "/Mesh/MixedElements";
  } else {
    tmp1 << stripFilename_(H5MeshFilename()) << ":/" << static_mesh_cycle_ << "/Mesh/MixedElements";
  }

  DataItem.addContent(tmp1.str());
  topo.addChild(DataItem);

  return topo;
}


Teuchos::XMLObject HDF5_MPI::addXdmfGeo_(const int cycle)
{
  std::stringstream tmp;
  std::stringstream tmp1;

  Teuchos::XMLObject geo("Geometry");
  geo.addAttribute("Name", "geo");
  geo.addAttribute("Type", "XYZ");

  Teuchos::XMLObject DataItem("DataItem");
  DataItem.addAttribute("DataType", "Float");
  tmp1 << NumNodes() << " " << " 3";
  DataItem.addAttribute("Dimensions", tmp1.str());
  DataItem.addAttribute("Format", "HDF");
  if (dynamic_mesh_) {
    tmp << stripFilename_(H5MeshFilename()) << ":/" << cycle << "/Mesh/Nodes";
  } else {
    tmp << stripFilename_(H5MeshFilename()) << ":/" << static_mesh_cycle_ << "/Mesh/Nodes";
  }
  DataItem.addContent(tmp.str());
  geo.addChild(DataItem);

  return geo;
}


void HDF5_MPI::writeXdmfVisitGrid_(std::string filename)
{
  // Create xmlObject grid
  Teuchos::XMLObject xi_include("xi:include");
  xi_include.addAttribute("href", stripFilename_(filename));
  xi_include.addAttribute("xpointer", "xpointer(//Xdmf/Domain/Grid)");

  // Step through xmlobject visit to find /domain/grid
  for (auto& xmf : xmlVisit_) { 
    Teuchos::XMLObject node;
    node = findGridNode_(xmf);
    node.addChild(xi_include);
  }
}


void HDF5_MPI::writeXdmfMeshVisitGrid_(std::string filename)
{
  // Create xmlObject grid
  Teuchos::XMLObject xi_include("xi:include");
  xi_include.addAttribute("href", stripFilename_(filename));
  xi_include.addAttribute("xpointer", "xpointer(//Xdmf/Domain/Grid)");

  // Step through xmlobject visit to find /domain/grid
  Teuchos::XMLObject node;
  node = findGridNode_(xmlMeshVisit_);

  // Add new grid to xmlobject visit
  node.addChild(xi_include);
}


Teuchos::XMLObject HDF5_MPI::findGridNode_(Teuchos::XMLObject xmlobject)
{
  Teuchos::XMLObject node, tmp;

  // Step down to child tag==Domain
  for (int i = 0; i < xmlobject.numChildren(); i++) {
    if (xmlobject.getChild(i).getTag() == "Domain") {
      node = xmlobject.getChild(i);
    }
  }

  // Step down to child tag==Grid and Attribute(GridType==Collection)
  for (int i = 0; i < node.numChildren(); i++) {
    tmp = node.getChild(i);
    if (tmp.getTag() == "Grid" && tmp.hasAttribute("GridType")) {
      if (tmp.getAttribute("GridType") == "Collection") {
        return tmp;
      }
    }
  }

  // TODO(barker): return some error indicator
  return node;
}


Teuchos::XMLObject HDF5_MPI::findMeshNode_(Teuchos::XMLObject xmlobject)
{
  Teuchos::XMLObject node, tmp;

  // Step down to child tag==Domain
  for (int i = 0; i < xmlobject.numChildren(); i++) {
    if (xmlobject.getChild(i).getTag() == "Domain") {
      node = xmlobject.getChild(i);
      break;
    }
  }

  // Step down to child tag==Grid and Attribute(Name==Mesh)
  for (int i = 0; i < node.numChildren(); i++) {
    tmp = node.getChild(i);
    if (tmp.getTag() == "Grid" && tmp.hasAttribute("Name")) {
      if (tmp.getAttribute("Name") == "Mesh") {
        return tmp;
      }
    }
  }

  // TODO(barker): return some error indicator
  return node;
}


Teuchos::XMLObject HDF5_MPI::addXdmfAttribute_(std::string varname,
                                               std::string location,
                                               int length,
                                               std::string h5path)
{
  Teuchos::XMLObject attribute("Attribute");
  attribute.addAttribute("Name", varname);
  attribute.addAttribute("Type", "Scalar");
  attribute.addAttribute("Center", location);

  Teuchos::XMLObject DataItem("DataItem");
  DataItem.addAttribute("Format", "HDF");
  DataItem.addInt("Dimensions", length);
  DataItem.addAttribute("DataType", "Float");
  std::stringstream tmp;
  tmp << stripFilename_(H5DataFilename()) << ":" << h5path;
  DataItem.addContent(tmp.str());
  attribute.addChild(DataItem);

  return attribute;
}


std::string HDF5_MPI::stripFilename_(std::string filename)
{
  std::stringstream ss(filename);
  std::string name;
  // strip for linux/unix/mac directory names
  char delim('/');
  while(std::getline(ss, name, delim)) {};
  // strip for windows directory names
  // delim='\\';
  // while(std::getline(ss, name, delim)) {}

  return name;
}

std::string HDF5_MPI::xdmfHeader_ = "<?xml version=\"1.0\" ?>\n<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";

}  // namespace Amanzi
