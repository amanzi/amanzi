/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

#include <fstream>
#include "HDF5Mesh.hh"

namespace Amanzi {

void
HDF5::createMeshFile(AmanziMesh::Mesh& mesh_maps, std::string filename)
{
  hid_t file, group, dataspace, dataset;
  hsize_t dimsf[2];
  std::string h5Filename, xmfFilename;

  // build h5 filename
  h5Filename = filename;

  file = H5Fcreate(h5Filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  if (file < 0) {
    Errors::Message message("HDF5:: error creating mesh file");
    Exceptions::amanzi_throw(message);
  }

  // TODO(barker): should have some more error handling here!
  group = H5Gcreate(file, "/Mesh", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  // get num_nodes, num_cells
  int num_nodes = mesh_maps.getNumEntities(AmanziMesh::Entity_kind::NODE, AmanziMesh::Parallel_kind::OWNED);
  int num_elems = mesh_maps.getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);

  // get coords
  double* nodes = new double[num_nodes * 3];

  AmanziGeometry::Point xc;
  for (int i = 0; i < num_nodes; i++) {
    xc = mesh_maps.getNodeCoordinate(i);
    nodes[i * 3] = xc[0];
    nodes[i * 3 + 1] = xc[1];
    nodes[i * 3 + 2] = xc[2];
  }

  // write out coords
  // TODO(barker): add error handling: can't create/write
  dimsf[0] = num_nodes;
  dimsf[1] = 3;
  dataspace = H5Screate_simple(2, dimsf, NULL);
  dataset =
    H5Dcreate(group, "Nodes", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, nodes);
  delete[] nodes;
  H5Dclose(dataset);
  H5Sclose(dataspace);

  // get number of nodes per element and element type
  ctype_ = mesh_maps.getCellType(0);
  cname_ = to_string(ctype_);
  unsigned int cellid = 0;
  auto nodeids = mesh_maps.getCellNodes(cellid);
  conn_ = nodeids.size();

  // get connectivity
  int* ielem = new int[num_elems * conn_];
  //std::vector<unsigned int> xh(8);

  for (unsigned int i = 0; i < num_elems; i++) {
    //mesh_maps.cell_to_nodes(i, xh.begin(), xh.end());
    auto nodeids = mesh_maps.getCellNodes(i);
    for (int j = 0; j < conn_; j++) { ielem[i * conn_ + j] = nodeids[j]; }
  }

  // write out connectivity
  // TODO(barker): add error handling: can't create/write
  dimsf[0] = num_elems;
  dimsf[1] = conn_;
  dataspace = H5Screate_simple(2, dimsf, NULL);
  dataset =
    H5Dcreate(group, "Elements", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, ielem);
  delete[] ielem;
  H5Dclose(dataset);
  H5Sclose(dataspace);

  // close file
  H5Fclose(file);

  // Store information
  setH5MeshFilename(h5Filename);
  setNumNodes(num_nodes);
  setNumElems(num_elems);

  // Create and write out accompanying Xdmf file
  if (TrackXdmf()) {
    xmfFilename = filename + ".xmf";
    createXdmfMesh_(xmfFilename);
  }
}

void
HDF5::createDataFile(std::string soln_filename)
{
  std::string h5filename, PVfilename, Vfilename;
  hid_t file;

  // ?? input mesh filename or grab global mesh filename
  // ->assumes global name exists!!
  // build h5 filename
  h5filename = soln_filename;

  file = H5Fcreate(h5filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  if (file < 0) {
    Errors::Message message("HDF5:: error creating data file");
    Exceptions::amanzi_throw(message);
  }

  // close file
  H5Fclose(file);

  // Store filenames
  setH5DataFilename(h5filename);
  if (TrackXdmf()) {
    setxdmfParaviewFilename(soln_filename + "PV.xmf");
    setxdmfVisitFilename(soln_filename + "V.xmf");

    // start xmf files xmlObjects stored inside functions
    createXdmfParaview_();
    createXdmfVisit_();
  }
}

void
HDF5::createTimestep(const double time, const int iteration)
{
  std::ofstream of;

  if (TrackXdmf()) {
    // create single step xdmf file
    Teuchos::XMLObject tmp("Xdmf");
    tmp.addChild(addXdmfHeaderLocal_(time));
    std::stringstream filename;
    filename << H5DataFilename() << "." << iteration << ".xmf";
    of.open(filename.str().c_str());
    of << tmp << std::endl;
    of.close();
    setxdmfStepFilename(filename.str());

    // update ParaView and VisIt xdmf files
    // TODO(barker): how to get to grid collection node, rather than root???
    writeXdmfParaviewGrid_(filename.str(), time, iteration);
    writeXdmfVisitGrid_(filename.str());
    // TODO(barker): where to write out depends on where the root node is
    // ?? how to terminate stream or switch to new file out??
    of.open(xdmfParaviewFilename().c_str());
    of << xmlParaview();
    of.close();
    of.open(xdmfVisitFilename().c_str());
    of << xmlVisit();
    of.close();

    // TODO(barker): where to add time & iteration information in h5 data file?
    // Store information
    xmlStep_ = tmp;
  }
  setIteration(iteration);
}

void
HDF5::endTimestep()
{
  if (TrackXdmf()) {
    std::ofstream of;
    std::stringstream filename;
    filename << H5DataFilename() << "." << Iteration() << ".xmf";
    of.open(filename.str().c_str());
    of << xmlStep();
    of.close();
  }
}

void
HDF5::writeCellData(const Epetra_Vector& x, const std::string varname)
{
  writeFieldData_(x, varname, "Cell");
}

void
HDF5::writeNodeData(const Epetra_Vector& x, const std::string varname)
{
  writeFieldData_(x, varname, "Node");
}

void
HDF5::writeData(const Epetra_Vector& x, const std::string varname)
{
  writeFieldData_(x, varname, "None");
}
void
HDF5::readData(Epetra_Vector& x, const std::string varname)
{
  readFieldData_(x, varname);
}


void
HDF5::createXdmfMesh_(const std::string xmfFilename)
{
  // TODO(barker): add error handling: can't open/write
  Teuchos::XMLObject mesh("Xdmf");

  // build xml object
  mesh.addChild(addXdmfHeaderLocal_(0.0));

  // write xmf
  std::ofstream of(xmfFilename.c_str());
  of << HDF5::xdmfHeader_ << mesh << std::endl;
  of.close();
}

void
HDF5::createXdmfParaview_()
{
  // TODO(barker): add error handling: can't open/write
  Teuchos::XMLObject xmf("Xdmf");
  xmf.addAttribute("xmlns:xi", "http://www.w3.org/2001/XInclude");
  xmf.addAttribute("Version", "2.0");

  // build xml object
  Teuchos::XMLObject header;
  header = addXdmfHeaderGlobal_();
  xmf.addChild(header);
  Teuchos::XMLObject node;
  node = findGridNode_(xmf);
  node.addChild(addXdmfTopo_());
  node.addChild(addXdmfGeo_());

  // write xmf
  std::ofstream of(xdmfParaviewFilename().c_str());
  of << HDF5::xdmfHeader_ << xmf << std::endl;
  of.close();

  // Store ParaView XMLObject
  xmlParaview_ = xmf;
}

void
HDF5::createXdmfVisit_()
{
  // TODO(barker): add error handling: can't open/write
  Teuchos::XMLObject xmf("Xdmf");
  xmf.addAttribute("xmlns:xi", "http://www.w3.org/2001/XInclude");
  xmf.addAttribute("Version", "2.0");

  // build xml object
  xmf.addChild(addXdmfHeaderGlobal_());

  // write xmf
  std::ofstream of(xdmfVisitFilename().c_str());
  of << HDF5::xdmfHeader_ << xmf << std::endl;
  of.close();

  // Store VisIt XMLObject
  xmlVisit_ = xmf;
}

Teuchos::XMLObject
HDF5::addXdmfHeaderGlobal_()
{
  Teuchos::XMLObject domain("Domain");

  Teuchos::XMLObject grid("Grid");
  grid.addAttribute("GridType", "Collection");
  grid.addAttribute("CollectionType", "Temporal");
  domain.addChild(grid);

  return domain;
}

Teuchos::XMLObject
HDF5::addXdmfHeaderLocal_(const double value)
{
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

Teuchos::XMLObject
HDF5::addXdmfTopo_()
{
  std::stringstream tmp, tmp1;

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

  topo.addChild(DataItem);

  return topo;
}

Teuchos::XMLObject
HDF5::addXdmfGeo_()
{
  std::string tmp;
  std::stringstream tmp1;

  Teuchos::XMLObject geo("Geometry");
  geo.addAttribute("Name", "geo");
  geo.addAttribute("Type", "XYZ");

  Teuchos::XMLObject DataItem("DataItem");
  DataItem.addAttribute("DataType", "Float");
  tmp1 << NumNodes() << " 3";
  DataItem.addAttribute("Dimensions", tmp1.str());
  DataItem.addAttribute("Format", "HDF");
  tmp = H5MeshFilename() + ":/Mesh/Nodes";
  DataItem.addContent(tmp);
  geo.addChild(DataItem);

  return geo;
}

void
HDF5::writeXdmfParaviewGrid_(std::string filename, const double time, const int iteration)
{
  // Create xmlObject grid

  Teuchos::XMLObject grid("Grid");
  grid.addInt("Name", iteration);
  grid.addAttribute("GridType", "Uniform");

  // TODO(barker): update topo/geo name if mesh is evolving
  Teuchos::XMLObject topo("Topology");
  topo.addAttribute("Reference", "//Topology[@Name='topo']");
  grid.addChild(topo);
  Teuchos::XMLObject geo("Geometry");
  geo.addAttribute("Reference", "//Geometry[@Name='geo']");
  grid.addChild(geo);

  Teuchos::XMLObject timeValue("Time");
  timeValue.addDouble("Value", time);
  grid.addChild(timeValue);

  Teuchos::XMLObject xi_include("xi:include");
  xi_include.addAttribute("href", filename);
  xi_include.addAttribute("xpointer", "xpointer(//Xdmf/Domain/Grid/Attribute)");
  grid.addChild(xi_include);

  // Step through paraview xmlobject to find /domain/grid
  Teuchos::XMLObject node;
  node = findGridNode_(xmlParaview_);

  // Add new grid to xmlobject paraview
  node.addChild(grid);
}

void
HDF5::writeXdmfVisitGrid_(std::string filename)
{
  // Create xmlObject grid
  Teuchos::XMLObject xi_include("xi:include");
  xi_include.addAttribute("href", filename);
  xi_include.addAttribute("xpointer", "xpointer(//Xdmf/Domain/Grid)");

  // Step through xmlobject visit to find /domain/grid
  Teuchos::XMLObject node;
  node = findGridNode_(xmlVisit_);

  // Add new grid to xmlobject visit
  node.addChild(xi_include);
}

Teuchos::XMLObject
HDF5::findGridNode_(Teuchos::XMLObject xmlobject)
{
  Teuchos::XMLObject node, tmp;

  // Step down to child tag==Domain
  for (int i = 0; i < xmlobject.numChildren(); i++) {
    if (xmlobject.getChild(i).getTag() == "Domain") { node = xmlobject.getChild(i); }
  }

  // Step down to child tag==Grid and Attribute(GridType==Collection)
  for (int i = 0; i < node.numChildren(); i++) {
    tmp = node.getChild(i);
    if (tmp.getTag() == "Grid" && tmp.hasAttribute("GridType")) {
      if (tmp.getAttribute("GridType") == "Collection") { return tmp; }
    }
  }

  // TODO(barker): return some error indicator
  return node;
}

Teuchos::XMLObject
HDF5::findMeshNode_(Teuchos::XMLObject xmlobject)
{
  Teuchos::XMLObject node, tmp;

  // Step down to child tag==Domain
  for (int i = 0; i < xmlobject.numChildren(); i++) {
    if (xmlobject.getChild(i).getTag() == "Domain") { node = xmlobject.getChild(i); }
  }

  // Step down to child tag==Grid and Attribute(Name==Mesh)
  for (int i = 0; i < node.numChildren(); i++) {
    tmp = node.getChild(i);
    if (tmp.getTag() == "Grid" && tmp.hasAttribute("Name")) {
      if (tmp.getAttribute("Name") == "Mesh") { return tmp; }
    }
  }

  // TODO(barker): return some error indicator
  return node;
}

void
HDF5::writeFieldData_(const Epetra_Vector& x, std::string varname, std::string loc)
{
  // write field data
  double* data;
  x.ExtractView(&data);
  int length = x.GlobalLength();
  hid_t file, group, dataspace, dataset;
  herr_t status;

  // TODO(barker): how to build path name?? probably still need iteration number
  std::stringstream h5path;
  h5path << "/" << varname;

  // TODO(barker): add error handling: can't write/create
  file = H5Fopen(H5DataFilename().c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
  if (file < 0) {
    Errors::Message message("HDF5::writeFieldData_ error opening data file to write field data");
    Exceptions::amanzi_throw(message);
  }

  if (TrackXdmf()) {
    // Check if varname group exists; if not, create it
    htri_t exists = H5Lexists(file, h5path.str().c_str(), H5P_DEFAULT);
    if (exists) {
      std::cout << "  WRITE>> opening group:" << h5path.str() << std::endl;
      group = H5Gopen(file, h5path.str().c_str(), H5P_DEFAULT);
    } else {
      std::cout << "  WRITE>> creating group:" << h5path.str() << std::endl;
      group = H5Gcreate(file, h5path.str().c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }
    h5path << "/" << Iteration();
    status = H5Gclose(group);
  }

  hsize_t hs = static_cast<hsize_t>(length);
  dataspace = H5Screate_simple(1, &hs, NULL);
  dataset = H5Dcreate(file,
                      h5path.str().c_str(),
                      H5T_NATIVE_DOUBLE,
                      dataspace,
                      H5P_DEFAULT,
                      H5P_DEFAULT,
                      H5P_DEFAULT);
  status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
  if (status < 0) {
    Errors::Message message("HDF5:: error writing field to data file");
    Exceptions::amanzi_throw(message);
  }
  status = H5Dclose(dataset);
  status = H5Sclose(dataspace);
  status = H5Fclose(file);

  // TODO(barker): add error handling: can't write
  if (TrackXdmf() and loc != "None") {
    // TODO(barker): get grid node, node.addChild(addXdmfAttribute)
    Teuchos::XMLObject node = findMeshNode_(xmlStep());
    node.addChild(addXdmfAttribute_(varname, loc, length, h5path.str()));
  }
}

void
HDF5::readFieldData_(Epetra_Vector& x, std::string varname)
{
  char* h5path = new char[varname.size() + 1];
  strcpy(h5path, varname.c_str());
  hid_t file = H5Fopen(H5DataFilename().c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
  int localdims[2];
  localdims[0] = x.MyLength();
  localdims[1] = 1;
  std::vector<int> myidx(localdims[0], 0);
  int start = 0;
  for (int i = 0; i < localdims[0]; i++) myidx[i] = i + start;
  double* data;
  x.ExtractView(&data);

  hid_t dataset = H5Dopen(file, h5path, H5P_DEFAULT);
  H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
  x.ReplaceMyValues(localdims[0], &data[0], &myidx[0]);

  H5Dclose(dataset);
  H5Fclose(file);
}

Teuchos::XMLObject
HDF5::addXdmfAttribute_(std::string varname, std::string location, int length, std::string h5path)
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
  tmp << H5DataFilename() << ":" << h5path;
  DataItem.addContent(tmp.str());
  attribute.addChild(DataItem);

  return attribute;
}


std::string HDF5::xdmfHeader_ =
  "<?xml version=\"1.0\" ?>\n<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";

} // namespace Amanzi
