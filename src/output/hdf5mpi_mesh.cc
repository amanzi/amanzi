#include "hdf5mpi_mesh.hh"

//TODO(barker): clean up debugging output
//TODO(barker): check that close file is always getting called
//TODO(barker): add error handling where appropriate
//TODO(barker): clean up formating

namespace Amanzi {
  

HDF5_MPI::HDF5_MPI(const Epetra_MpiComm &comm)
: viz_comm_(comm)
{
  viz_comm_ = comm;
  info_ = MPI_INFO_NULL;
  IOconfig_.numIOgroups = 1;
  IOconfig_.commIncoming = comm.Comm();
  parallelIO_IOgroup_init(&IOconfig_, &IOgroup_);
}

HDF5_MPI::HDF5_MPI(const Epetra_MpiComm &comm, std::string dataFilename)
: viz_comm_(comm), H5DataFilename_(dataFilename)
{
  viz_comm_ = comm;
  H5DataFilename_ = dataFilename; 
  info_ = MPI_INFO_NULL;
  IOconfig_.numIOgroups = 1;
  IOconfig_.commIncoming = comm.Comm();
  parallelIO_IOgroup_init(&IOconfig_, &IOgroup_);
}

HDF5_MPI::~HDF5_MPI()
{
  parallelIO_IOgroup_cleanup(&IOgroup_);
}

void HDF5_MPI::createMeshFile(const AmanziMesh::Mesh &mesh_maps, std::string filename)
{

  hid_t file, group, dataspace, dataset;
  herr_t status;
  hsize_t dimsf[2];
  std::string h5Filename, xmfFilename;
  int globaldims[2], localdims[2];
  int *ids;

  // build h5 filename
  h5Filename = filename;
  //TODO(barker): verify with Markus whether he'll include h5 or I'll add it
  h5Filename.append(std::string(".h5"));

  // new parallel
  file = parallelIO_open_file(h5Filename.c_str(), &IOgroup_, FILE_CREATE);
  if (file < 0) {
    Errors::Message message("HDF5_MPI::createMeshFile - error creating mesh file");
    Exceptions::amanzi_throw(message);
  }
  
  // get num_nodes, num_cells
  const Epetra_Map &nmap = mesh_maps.node_map(false);
  int nnodes_local = nmap.NumMyElements();
  int nnodes_global = nmap.NumGlobalElements();
  const Epetra_Map &ngmap = mesh_maps.node_map(true);
  
  const Epetra_Map &cmap = mesh_maps.cell_map(false);
  int ncells_local = cmap.NumMyElements();
  int ncells_global = cmap.NumGlobalElements();
  
  // get space dimension
  int space_dim = mesh_maps.space_dimension();
  //AmanziGeometry::Point xc;
  //mesh_maps.node_get_coordinates(0, &xc);
  //unsigned int space_dim = xc.dim();

  // get coords
  double *nodes = new double[nnodes_local*3];
  globaldims[0] = nnodes_global;
  globaldims[1] = 3;
  localdims[0] = nnodes_local;
  localdims[1] = 3;

  std::vector<double> xc(space_dim);
  for (int i = 0; i < nnodes_local; i++) {
    mesh_maps.node_to_coordinates(i, xc.begin(), xc.end());
    // VisIt and ParaView require all mesh entities to be in 3D space
    nodes[i*3+0] = xc[0];
    nodes[i*3+1] = xc[1];
    if (space_dim == 3) {
      nodes[i*3+2] = xc[2];
    } else {
      nodes[i*3+2] = 0.0;
    }
  }

  // write out coords
  // TODO(barker): add error handling: can't create/write
  parallelIO_write_dataset(nodes, PIO_DOUBLE, 2, globaldims, localdims, file,
                           (char*)"Mesh/Nodes", &IOgroup_, 
                           NONUNIFORM_CONTIGUOUS_WRITE);
  delete nodes;
  
  // write out node map
  ids = new int[nmap.NumMyElements()];
  for (int i=0; i<nnodes_local; i++) {
    ids[i] = nmap.GID(i);
  }
  globaldims[1] = 1;
  localdims[1] = 1;
  parallelIO_write_dataset(ids, PIO_INTEGER, 2, globaldims, localdims, file,
                           (char*)"Mesh/NodeMap", &IOgroup_,
                           NONUNIFORM_CONTIGUOUS_WRITE);
  
  
  // get connectivity
  // nodes are written to h5 out of order, need info to map id to order in output
  int nnodes(nnodes_local);
  std::vector<int> nnodesAll(viz_comm_.NumProc(),0);
  viz_comm_.GatherAll(&nnodes, &nnodesAll[0], 1);
  int start(0);
  std::vector<int> startAll(viz_comm_.NumProc(),0);
  for (int i = 0; i < viz_comm_.MyPID(); i++){
    start += nnodesAll[i];
  }
  viz_comm_.GatherAll(&start, &startAll[0],1);
  
  AmanziMesh::Entity_ID_List nodeids;
  unsigned int cellid = 0;
  mesh_maps.cell_get_nodes(cellid,&nodeids);
  conn_ = nodeids.size();
  
  std::vector<int> gid(nnodes_global);
  std::vector<int> pid(nnodes_global);
  std::vector<int> lid(nnodes_global);
  for (int i=0; i<nnodes_global; i++) {
    gid[i] = ngmap.GID(i);
  }
  nmap.RemoteIDList(nnodes_global, &gid[0], &pid[0], &lid[0]);

  // determine size of connectivity vector 
  // element conn vector: elem_typeID elem_conn1 ... elem_connN
  // conn vector length = size_conn + 1
  // if polygon: elem_typeID num_nodes elem_conn1 ... elem_connN
  //             conn vector length = size_conn + 1 + 1
  // TODO(barker): make a list of cell types found, 
  //               if all the same then write out as a uniform mesh of that type
  int local_conn(0);
  AmanziMesh::Cell_type type;
  std::vector<int> each_conn(ncells_local);
  for (int i=0; i<ncells_local; i++) {
    mesh_maps.cell_get_nodes(i,&nodeids);
    each_conn[i] = nodeids.size();
    local_conn += each_conn[i]+1;  // add 1 for elem_typeID
    type = mesh_maps.cell_get_type(i); 
    if ( getCellTypeID_(type) == 3) local_conn += 1; // add 1 if polygon
  }
  std::vector<int> local_connAll(viz_comm_.NumProc(),0);
  viz_comm_.GatherAll(&local_conn, &local_connAll[0], 1);
  int total_conn(0);
  for (int i=0; i<viz_comm_.NumProc(); i++) {
    total_conn += local_connAll[i];
  }
  
  int *cells = new int[local_conn];
  globaldims[0] = total_conn; 
  globaldims[1] = 1;
  localdims[0] = local_conn;
  localdims[1] = 1;
  
  // get local element connectivities
  // nodeIDs need to be mapped to output IDs
  int idx = 0;
  for (int i=0; i<ncells_local; i++) {
    int conn_len(each_conn[i]);
    std::vector<unsigned int> xe(conn_len);
    mesh_maps.cell_to_nodes(i, xe.begin(), xe.end());
    // store cell type id
    type = mesh_maps.cell_get_type(i);
    cells[idx] = getCellTypeID_(type); 
    idx++;
    // TODO(barker): this shouldn't be a hardcoded value
    if (type == 3) {
      cells[idx+1] = each_conn[i];
      idx++;
    } 
    // store mapped node ids for connectivity
    for (int j = 0; j < conn_len; j++) {
      if (nmap.MyLID(xe[j])) {
        cells[idx+j] = xe[j] + startAll[viz_comm_.MyPID()];
      } else {
	cells[idx+j] = lid[xe[j]] + startAll[pid[xe[j]]];
      }
    }
    idx += conn_len;
  }
  // write out connectivity
  parallelIO_write_dataset(cells, PIO_INTEGER, 2, globaldims, localdims, file, 
                           (char*)"Mesh/MixedElements", &IOgroup_,
                           NONUNIFORM_CONTIGUOUS_WRITE);
  delete cells;
  
  // write out cell map
  ids = new int[cmap.NumMyElements()];
  for (int i=0; i<ncells_local; i++) {
    ids[i] = cmap.GID(i);
  }

  globaldims[0] = ncells_global; 
  globaldims[1] = 1;
  localdims[0] = ncells_local;
  localdims[1] = 1;
  parallelIO_write_dataset(ids, PIO_INTEGER, 2, globaldims, localdims, file, 
                           (char*)"Mesh/ElementMap", &IOgroup_,
                           NONUNIFORM_CONTIGUOUS_WRITE);
  delete ids;
  
  // close file
  parallelIO_close_file(file, &IOgroup_);

  // Store information
  setH5MeshFilename(h5Filename);
  setNumNodes(nnodes_global);
  setNumElems(ncells_global);
  setConnLength(total_conn);
  //TODO(barker): store the connectivity length for mixed meshes, anything else?

  // Create and write out accompanying Xdmf file
  if (TrackXdmf() && viz_comm_.MyPID() == 0) {
    //TODO(barker): if implement type tracking, then update this as needed
    ctype_ = mesh_maps.cell_get_type(0); 
    cname_ = "Mixed"; //AmanziMesh::Data::type_to_name(ctype_);
    xmfFilename = filename + ".xmf";
    createXdmfMesh_(xmfFilename);
  }

}

void HDF5_MPI::createDataFile(std::string soln_filename) {

  std::string h5filename, PVfilename, Vfilename;
  hid_t file;
  herr_t status;

  // ?? input mesh filename or grab global mesh filename
  // ->assumes global name exists!!
  // build h5 filename
  h5filename = soln_filename;
  //TODO:barker - verify with Markus whether he'll include h5 or I'll add it
  h5filename.append(std::string(".h5"));
  
  // new parallel
  file = parallelIO_open_file(h5filename.c_str(), &IOgroup_, FILE_CREATE);
  if (file < 0) {
    Errors::Message message("HDF5_MPI::createDataFile - error creating data file");
    Exceptions::amanzi_throw(message);
  }

  // close file
  parallelIO_close_file(file, &IOgroup_);

  // Store filenames
  setH5DataFilename(h5filename);
  if (TrackXdmf() && viz_comm_.MyPID() == 0) {
    setxdmfVisitFilename(soln_filename + ".VisIt.xmf");
    // start xmf files xmlObjects stored inside functions
    createXdmfVisit_();
  }
}

void HDF5_MPI::createTimestep(const double time, const int iteration) {

  std::ofstream of;

  if (TrackXdmf() && viz_comm_.MyPID() == 0) {
    // create single step xdmf file
    Teuchos::XMLObject tmp("Xdmf");
    tmp.addChild(addXdmfHeaderLocal_(time));
    std::stringstream filename;
    filename << H5DataFilename() << "." << iteration << ".xmf";
    of.open(filename.str().c_str());
    of << tmp << endl;
    of.close();
    setxdmfStepFilename(filename.str());

    // update VisIt xdmf files
    // TODO(barker): how to get to grid collection node, rather than root???
    writeXdmfVisitGrid_(filename.str());
    // TODO(barker): where to write out depends on where the root node is
    // ?? how to terminate stream or switch to new file out??
    of.open(xdmfVisitFilename().c_str());
    of << xmlVisit();
    of.close();

    // TODO(barker): where to add time & iteration information in h5 data file?
    // Store information
    xmlStep_ = tmp;
  }
  setIteration(iteration);
  setTime(time);
}

void HDF5_MPI::endTimestep() {
  if (TrackXdmf() && viz_comm_.MyPID() == 0) {
    std::ofstream of;
    std::stringstream filename;
    filename << H5DataFilename() << "." << Iteration() << ".xmf";
    of.open(filename.str().c_str());
    of << xmlStep();
    of.close();
  }
}

void HDF5_MPI::writeAttrString(const std::string value, const std::string attrname)
{
  
  hid_t file, fid, dataspace_id, attribute_id, atype;
  hsize_t dims[1]={1};
  herr_t status;
  
  /* Open the file */
  hid_t acc_tpl1 = H5Pcreate (H5P_FILE_ACCESS);
  herr_t ret = H5Pset_fapl_mpio(acc_tpl1, viz_comm_.Comm(), info_);
  file = H5Fopen(H5DataFilename_.c_str(), H5F_ACC_RDWR, acc_tpl1);
    
  if (file < 0) {
    Errors::Message message("HDF5_MPI::writeAttrString - error opening data file to write attribute");
    Exceptions::amanzi_throw(message);
  }
  
  /* Create a dataset attribute. */
  dataspace_id = H5Screate(H5S_SCALAR);
  atype = H5Tcopy(H5T_C_S1);
  H5Tset_size(atype, strlen(value.c_str())+1);
  H5Tset_strpad(atype,H5T_STR_NULLTERM);
  H5Tset_order(atype,H5T_ORDER_NONE);
  attribute_id = H5Acreate (file, attrname.c_str(), atype, 
                            dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
  
  /* Write the attribute data. */
  status = H5Awrite(attribute_id, atype, value.c_str());
  
  /* Close everything. */
  status = H5Aclose(attribute_id);  
  status = H5Sclose(dataspace_id);  
  status = H5Tclose(atype);  
  status = H5Pclose(acc_tpl1);
  status = H5Fclose(file);
}
  
  
void HDF5_MPI::writeAttrReal(double value, const std::string attrname)
{
  
  hid_t file, fid, dataspace_id, attribute_id;
  hsize_t dims;
  herr_t status;
  
  /* Setup dataspace */
  dims = 1;
  dataspace_id = H5Screate_simple(1, &dims, NULL);
  
  /* Open the file */
  hid_t acc_tpl1 = H5Pcreate (H5P_FILE_ACCESS);
  herr_t ret = H5Pset_fapl_mpio(acc_tpl1, viz_comm_.Comm(), info_);
  file = H5Fopen(H5DataFilename_.c_str(), H5F_ACC_RDWR, acc_tpl1);
    
  if (file < 0) {
    Errors::Message message("HDF5_MPI::writeAttrReal - error opening data file to write attribute");
    Exceptions::amanzi_throw(message);
  }
    
  /* Create a dataset attribute. */
  attribute_id = H5Acreate (file, attrname.c_str(), H5T_NATIVE_DOUBLE, 
                            dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
  
  /* Write the attribute data. */
  status = H5Awrite(attribute_id, H5T_NATIVE_DOUBLE, &value);
  
  /* Close everything. */
  status = H5Aclose(attribute_id);  
  status = H5Sclose(dataspace_id);  
  status = H5Pclose(acc_tpl1);
  status = H5Fclose(file);
}

  void HDF5_MPI::writeAttrReal(double value, const std::string attrname, std::string h5path)
  {
    
    hid_t file, dataset, fid, dataspace_id, attribute_id;
    hsize_t dims;
    herr_t status;
    char *DSpath = new char [h5path.size()+1];
    
    strcpy(DSpath,h5path.c_str());
    /* Setup dataspace */
    dims = 1;
    dataspace_id = H5Screate_simple(1, &dims, NULL);
    
    /* Open the file */
    hid_t acc_tpl1 = H5Pcreate (H5P_FILE_ACCESS);
    herr_t ret = H5Pset_fapl_mpio(acc_tpl1, viz_comm_.Comm(), info_);
    file = H5Fopen(H5DataFilename_.c_str(), H5F_ACC_RDWR, acc_tpl1);
    if (file < 0) {
      Errors::Message message("HDF5_MPI::writeAttrReal2 - error opening data file to write attribute");
      Exceptions::amanzi_throw(message);
    }
    
    dataset = H5Dopen(file, DSpath, H5P_DEFAULT);
    if (dataset < 0) {
      Errors::Message message("HDF5_MPI::writeAttrReal2 - error opening dataset in file to write attribute");
      Exceptions::amanzi_throw(message);
    }
    
    /* Create a dataset attribute. */
    attribute_id = H5Acreate (dataset, attrname.c_str(), H5T_NATIVE_DOUBLE, 
                              dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
    
    /* Write the attribute data. */
    status = H5Awrite(attribute_id, H5T_NATIVE_DOUBLE, &value);
    
    /* Close everything. */
    status = H5Dclose(dataset);
    status = H5Aclose(attribute_id);  
    status = H5Sclose(dataspace_id);  
    status = H5Pclose(acc_tpl1);
    status = H5Fclose(file);
  }
  
void HDF5_MPI::writeAttrInt(int value, const std::string attrname)
{
  
  hid_t file, dataspace_id, attribute_id;
  hsize_t dims;
  herr_t status;

  /* Setup dataspace */
  dims = 1;
  dataspace_id = H5Screate_simple(1, &dims, NULL);
  
  /* Open the file */
  //file = parallelIO_open_file(H5DataFilename_.c_str(), &IOgroup_,
  //                            FILE_READWRITE);
  //iofile_t *curr_file = IOgroup_.file[file];
  hid_t acc_tpl1 = H5Pcreate (H5P_FILE_ACCESS);
  herr_t ret = H5Pset_fapl_mpio(acc_tpl1, viz_comm_.Comm(), info_);
  file = H5Fopen(H5DataFilename_.c_str(), H5F_ACC_RDWR, acc_tpl1);
    
  if (file < 0) {
    Errors::Message message("HDF5_MPI::writeAttrInt - error opening data file to write attribute");
    Exceptions::amanzi_throw(message);
  }
    
  /* Create/Write a dataset attribute. */
  attribute_id = H5Acreate (file, attrname.c_str(), H5T_NATIVE_INT,
                            dataspace_id, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Awrite(attribute_id, H5T_NATIVE_INT, &value);
  
  /* Set the attribute data. */
  //H5LTset_attribute_int(curr_file->fid, "/", attrname.c_str(), &value, 1);
  //H5LTset_attribute_int(file, "/", attrname.c_str(), &value, 1);
  
  /* Close everything. */
  status = H5Aclose(attribute_id);  
  status = H5Sclose(dataspace_id); 
  status = H5Pclose(acc_tpl1);
  status = H5Fclose(file);
  //parallelIO_close_file(file, &IOgroup_);
  
}

void HDF5_MPI::readAttrString(std::string &value, const std::string attrname)
{
  
  hid_t file, fid, attribute_id, atype, root;
  herr_t status;
  
  /* Open the file */
  hid_t acc_tpl1 = H5Pcreate (H5P_FILE_ACCESS);
  herr_t ret = H5Pset_fapl_mpio(acc_tpl1, viz_comm_.Comm(), info_);
  file = H5Fopen(H5DataFilename_.c_str(), H5F_ACC_RDWR, acc_tpl1);
  root = H5Gopen(file, "/", H5P_DEFAULT);
    
  if (file < 0) {
    Errors::Message message("HDF5_MPI::readAttrString - error opening data file to write attribute");
    Exceptions::amanzi_throw(message);
  }
  
  atype = H5Tcopy(H5T_C_S1);
  attribute_id = H5Aopen(root, attrname.c_str(), H5P_DEFAULT);
  hid_t ftype = H5Aget_type(attribute_id);
  size_t size = H5Tget_size(ftype);
  H5Tset_size(atype, size);
  
  /* Write the attribute data. */
  char string_out[size];
  status = H5Aread(attribute_id, atype, string_out);
  value = string_out;

  /* Close everything. */ 
  status = H5Tclose(atype);    
  status = H5Tclose(ftype);  
  status = H5Pclose(acc_tpl1);
  status = H5Aclose(attribute_id); 
  status = H5Gclose(root);
  status = H5Fclose(file);
}
  
void HDF5_MPI::readAttrReal(double &value, const std::string attrname)
{
    
  hid_t file, fid, attribute_id;
  herr_t status;
    
  /* Open the file */
  hid_t acc_tpl1 = H5Pcreate (H5P_FILE_ACCESS);
  herr_t ret = H5Pset_fapl_mpio(acc_tpl1, viz_comm_.Comm(), info_);
  file = H5Fopen(H5DataFilename_.c_str(), H5F_ACC_RDWR, acc_tpl1);
    
  if (file < 0) {
    Errors::Message message("HDF5_MPI::readAttrReal - error opening data file to write attribute");
    Exceptions::amanzi_throw(message);
  }
    
  /* Create a dataset attribute. */
  attribute_id = H5Aopen(file, attrname.c_str(), H5P_DEFAULT);
    
  /* Write the attribute data. */
  status = H5Aread(attribute_id, H5T_NATIVE_DOUBLE, &value);
    
  /* Close everything. */
  status = H5Aclose(attribute_id);  
  status = H5Pclose(acc_tpl1);
  status = H5Fclose(file);
  
}
  
void HDF5_MPI::readAttrInt(int &value, const std::string attrname)
  {
    
    hid_t file, fid, attribute_id;
    herr_t status;
    
    /* Open the file */
    hid_t acc_tpl1 = H5Pcreate (H5P_FILE_ACCESS);
    herr_t ret = H5Pset_fapl_mpio(acc_tpl1, viz_comm_.Comm(), info_);
    file = H5Fopen(H5DataFilename_.c_str(), H5F_ACC_RDWR, acc_tpl1);
    
    if (file < 0) {
      Errors::Message message("HDF5_MPI::readAttrInt - error opening data file to write attribute");
      Exceptions::amanzi_throw(message);
    }
    
    /* Create a dataset attribute. */
    attribute_id = H5Aopen(file, attrname.c_str(), H5P_DEFAULT);
    
    /* Write the attribute data. */
    status = H5Aread(attribute_id, H5T_NATIVE_INT, &value);
    
    /* Close everything. */
    status = H5Aclose(attribute_id);  
    status = H5Pclose(acc_tpl1);
    status = H5Fclose(file);
  }
  
void HDF5_MPI::writeDataReal(const Epetra_Vector &x, const std::string varname) 
{
  writeFieldData_(x, varname, PIO_DOUBLE, "NONE");
}

void HDF5_MPI::writeDataInt(const Epetra_Vector &x, const std::string varname)
{
  writeFieldData_(x, varname, PIO_INTEGER, "NONE");
}
  
void HDF5_MPI::writeCellDataReal(const Epetra_Vector &x,
                                 const std::string varname)
{
  writeFieldData_(x, varname, PIO_DOUBLE, "Cell");
}

void HDF5_MPI::writeCellDataInt(const Epetra_Vector &x,
                                const std::string varname)
{
  writeFieldData_(x, varname, PIO_INTEGER, "Cell");
}

void HDF5_MPI::writeNodeDataReal(const Epetra_Vector &x, 
                                 const std::string varname)
{
  writeFieldData_(x, varname, PIO_DOUBLE, "Node");
}

void HDF5_MPI::writeNodeDataInt(const Epetra_Vector &x, 
                                const std::string varname)
{
  writeFieldData_(x, varname, PIO_INTEGER, "Node");
}

  
void HDF5_MPI::writeFieldData_(const Epetra_Vector &x, std::string varname,
                               datatype_t type, std::string loc) {
  // write field data
  double *data;
  int err = x.ExtractView(&data);
  hid_t file, group, dataspace, dataset;
  herr_t status;

  int globaldims[2], localdims[2];
  globaldims[0] = x.GlobalLength();
  globaldims[1] = 1;
  localdims[0] = x.MyLength();
  localdims[1] = 1;

  // TODO(barker): how to build path name?? probably still need iteration number
  std::stringstream h5path;
  h5path << varname;

  // TODO(barker): add error handling: can't write/create

  file = parallelIO_open_file(H5DataFilename_.c_str(), &IOgroup_,
                              FILE_READWRITE);
  if (file < 0) {
    Errors::Message message("HDF5_MPI::writeFieldData_ - error opening data file to write field data");
    Exceptions::amanzi_throw(message);
  }

  if (TrackXdmf()){
    h5path << "/" << Iteration();
  }
  
  char *tmp;
  tmp = new char [h5path.str().size()+1];
  strcpy(tmp,h5path.str().c_str());
  parallelIO_write_dataset(data, type, 2, globaldims, localdims, file, tmp,
                           &IOgroup_, NONUNIFORM_CONTIGUOUS_WRITE);    
  parallelIO_close_file(file, &IOgroup_);

  // TODO(barker): add error handling: can't write
  if (TrackXdmf() ) {
    // write the time value as an attribute to this dataset
    writeAttrReal(Time(), "Time", h5path.str());
    if (viz_comm_.MyPID() == 0) {
      // TODO(barker): get grid node, node.addChild(addXdmfAttribute)
      Teuchos::XMLObject node = findMeshNode_(xmlStep());
      node.addChild(addXdmfAttribute_(varname, loc, globaldims[0], h5path.str()));
    }
  }
}
  
void HDF5_MPI::readData(Epetra_Vector &x, const std::string varname)
{
  readFieldData_(x, varname, PIO_DOUBLE);
}
  
void HDF5_MPI::readFieldData_(Epetra_Vector &x, std::string varname,
                               datatype_t type) {
  
  char *h5path = new char [varname.size()+1];
  strcpy(h5path,varname.c_str());
  int ndims;
  hid_t file = parallelIO_open_file(H5DataFilename_.c_str(), &IOgroup_,
                                    FILE_READONLY);
  if (file < 0) {
    Errors::Message message("HDF5_MPI::readFieldData_ - error opening data file to write field data");
    Exceptions::amanzi_throw(message);
  }
  parallelIO_get_dataset_ndims(&ndims, file, h5path, &IOgroup_);
  int  globaldims[ndims], localdims[ndims];
  parallelIO_get_dataset_dims(globaldims, file, h5path, &IOgroup_);
  localdims[0] = x.MyLength();
  localdims[1] = globaldims[1];
  std::vector<int> myidx(localdims[0],0);
  int start = 0;
  for (int i=0; i<localdims[0]; i++) myidx[i] = i+start;
  
  double *data = new double[localdims[0]*localdims[1]];
  parallelIO_read_dataset(data, type, ndims, globaldims, localdims,
                          file, h5path, &IOgroup_, NONUNIFORM_CONTIGUOUS_READ);
  x.ReplaceMyValues(localdims[0], &data[0], &myidx[0]);
  
  parallelIO_close_file(file, &IOgroup_); 
  delete data;
  
}
  
int HDF5_MPI::getCellTypeID_(AmanziMesh::Cell_type type) {
  
  //TODO(barker): how to return polyhedra?
  // cell type id's defined in Xdmf/include/XdmfTopology.h
  
  ASSERT (cell_valid_type(type));

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
      return 3; //for now same as polygon
    default:
      return 3; //unknown, for now same as polygon
    }

}
  
void HDF5_MPI::createXdmfMesh_(const std::string xmfFilename) {
  // TODO(barker): add error handling: can't open/write
  Teuchos::XMLObject mesh("Xdmf");

  // build xml object
  mesh.addChild(addXdmfHeaderLocal_(0.0));

  // write xmf
  std::ofstream of(xmfFilename.c_str());
  of << HDF5_MPI::xdmfHeader_ << mesh << endl;
  of.close();
}


void HDF5_MPI::createXdmfVisit_() {
  // TODO(barker): add error handling: can't open/write
  Teuchos::XMLObject xmf("Xdmf");
  xmf.addAttribute("xmlns:xi", "http://www.w3.org/2001/XInclude");
  xmf.addAttribute("Version", "2.0");

  // build xml object
  xmf.addChild(addXdmfHeaderGlobal_());

  // write xmf
  std::ofstream of(xdmfVisitFilename().c_str());
  of << HDF5_MPI::xdmfHeader_ << xmf << endl;
  of.close();

  // Store VisIt XMLObject
  xmlVisit_ = xmf;
}

Teuchos::XMLObject HDF5_MPI::addXdmfHeaderGlobal_() {

  Teuchos::XMLObject domain("Domain");

  Teuchos::XMLObject grid("Grid");
  grid.addAttribute("GridType", "Collection");
  grid.addAttribute("CollectionType", "Temporal");
  domain.addChild(grid);

  return domain;
}

Teuchos::XMLObject HDF5_MPI::addXdmfHeaderLocal_(const double value) {
  Teuchos::XMLObject domain("Domain");

  Teuchos::XMLObject grid("Grid");
  grid.addAttribute("Name", "Mesh");
  domain.addChild(grid);
  grid.addChild(addXdmfTopo_());
  grid.addChild(addXdmfGeo_());

  Teuchos::XMLObject time("Time");
  time.addDouble("Value", value);
  grid.addChild(time);

  return domain;
}

Teuchos::XMLObject HDF5_MPI::addXdmfTopo_() {
  std::stringstream tmp, tmp1;
  
  // TODO(barker): error checking if cname_ and conn_ haven't been checked
  // TODO(barker): error checking if cname_ is unknown
  /*
  Teuchos::XMLObject topo("Topology");
  topo.addAttribute("Type", cname_);
  topo.addInt("Dimensions", NumElems());
  topo.addAttribute("Name", "topo");

  Teuchos::XMLObject DataItem("DataItem");
  DataItem.addAttribute("DataType", "Int");
  tmp << NumElems() << " " << conn_;
  DataItem.addAttribute("Dimensions", tmp.str());
  DataItem.addAttribute("Format", "HDF");

  tmp1 << H5MeshFilename() << ":/Mesh/Elements";
  DataItem.addContent(tmp1.str());
   */

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
  
  tmp1 << stripFilename_(H5MeshFilename()) << ":/Mesh/MixedElements";
  DataItem.addContent(tmp1.str());
  topo.addChild(DataItem);

  return topo;
}

Teuchos::XMLObject HDF5_MPI::addXdmfGeo_() {
  std::string tmp;
  std::stringstream tmp1;

  Teuchos::XMLObject geo("Geometry");
  geo.addAttribute("Name", "geo");
  geo.addAttribute("Type", "XYZ");

  Teuchos::XMLObject DataItem("DataItem");
  DataItem.addAttribute("DataType", "Float");
  tmp1 << NumNodes() << " " << " 3";
  DataItem.addAttribute("Dimensions", tmp1.str());
  DataItem.addAttribute("Format", "HDF");
  tmp = stripFilename_(H5MeshFilename()) + ":/Mesh/Nodes";
  DataItem.addContent(tmp);
  geo.addChild(DataItem);

  return geo;
}

void HDF5_MPI::writeXdmfVisitGrid_(std::string filename) {

  // Create xmlObject grid
  Teuchos::XMLObject xi_include("xi:include");
  xi_include.addAttribute("href", stripFilename_(filename));
  xi_include.addAttribute("xpointer", "xpointer(//Xdmf/Domain/Grid)");

  // Step through xmlobject visit to find /domain/grid
  Teuchos::XMLObject node;
  node = findGridNode_(xmlVisit_);

  // Add new grid to xmlobject visit
  node.addChild(xi_include);
}

Teuchos::XMLObject HDF5_MPI::findGridNode_(Teuchos::XMLObject xmlobject) {

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

Teuchos::XMLObject HDF5_MPI::findMeshNode_(Teuchos::XMLObject xmlobject) {
  
  Teuchos::XMLObject node, tmp;
  
  // Step down to child tag==Domain
  for (int i = 0; i < xmlobject.numChildren(); i++) {
    if (xmlobject.getChild(i).getTag() == "Domain") {
      node = xmlobject.getChild(i);
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
                                           std::string h5path) {
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

std::string HDF5_MPI::stripFilename_(std::string filename)  {
  
  std::stringstream ss(filename);
  std::string name;
  // strip for linux/unix/mac directory names
  char delim('/');
  while(std::getline(ss, name, delim)) { }
  // strip for windows directory names
  //delim='\\';
  //while(std::getline(ss, name, delim)) { }

  return name;
}
  
std::string HDF5_MPI::xdmfHeader_ =
         "<?xml version=\"1.0\" ?>\n<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
  
} // close namespace Amanzi

