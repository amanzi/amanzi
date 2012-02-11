#ifndef HDF5MPI_MESH_HH_
#define HDF5MPI_MESH_HH_

#include <string>

//#ifdef HAVE_MOAB_MESH
//#include "Mesh_moab.hh"
//#endif
//#ifdef HAVE_STK_MESH
//#include "Mesh_STK.hh"
//#include "Mesh_maps_stk.hh"
//#endif

#include "Mesh.hh"
#include "Element_types.hh"
#include "errors.hh"
#include "Epetra_Vector.h"
#include "Teuchos_XMLObject.hpp"

extern "C" {
#include "hdf5.h"
#include "hdf5_hl.h"
#include "parallelIO.h"
};

namespace Amanzi
{

class HDF5_MPI {


 public:

  HDF5_MPI(const Epetra_MpiComm &comm);
  HDF5_MPI(const Epetra_MpiComm &comm, std::string dataFilename);
  ~HDF5_MPI(void);
  
  bool TrackXdmf() { return TrackXdmf_; }
  void setTrackXdmf(bool TrackXdmf) { TrackXdmf_ = TrackXdmf; }

  std::string H5MeshFilename() { return H5MeshFilename_; }
  void setH5MeshFilename(std::string H5MeshFilename) {
                         H5MeshFilename_ = H5MeshFilename;}
  std::string H5DataFilename() { return H5DataFilename_; }
  void setH5DataFilename(std::string H5DataFilename) {
                         H5DataFilename_ = H5DataFilename;}
  std::string xdmfVisitFilename() { return xdmfVisitFilename_; }
  void setxdmfVisitFilename(std::string xdmfVisitFilename) {
                            xdmfVisitFilename_ = xdmfVisitFilename;}
  std::string xdmfStepFilename() { return xdmfStepFilename_; }
  void setxdmfStepFilename(std::string xdmfStepFilename) {
                           xdmfStepFilename_ = xdmfStepFilename;}
  
  int NumNodes() { return NumNodes_;}
  void setNumNodes(int NumNodes) {NumNodes_ = NumNodes;}
  int NumElems() { return NumElems_;}
  void setNumElems(int NumElems) {NumElems_ = NumElems;}
  int ConnLength() { return ConnLength_;}
  void setConnLength(int ConnLength) {ConnLength_ = ConnLength;}
  int Iteration() { return Iteration_;}
  void setIteration(int Iteration) {Iteration_ = Iteration;}

  Teuchos::XMLObject xmlVisit() { return xmlVisit_; }
  Teuchos::XMLObject xmlStep() { return xmlStep_; }

  // Output mesh data to filename.h5 and filename.xmf
  void createMeshFile(const AmanziMesh::Mesh &mesh_Maps, std::string filename);
  //void createMeshFile(Mesh_maps_base &mesh_Maps, std::string filename);

  // Create h5 file for data output, create accompanying Xdmf files for Visit
  void createDataFile(std::string data_filename);
  // Adds time step attributes to VisIt Xdmf files.  Creates
  // individual Xdmf for the current step.
  // TODO(barker): The individual step file can be remove after VisIt updates.
  // TODO(barker): Consolidate into a singel Xdmf file, after VisIt updates.
  void createTimestep(const double time, const int iteration);
  void endTimestep();

  // Write attribute to HDF5 data file.
  void writeAttrReal(double value, const std::string attrname);
  void writeAttrInt(int value, const std::string attrname);
  void writeAttrString(const std::string value, const std::string attrname);
  void readAttrReal(double &value, const std::string attrname);
  void readAttrInt(int &value, const std::string attrname);
  void readAttrString(std::string &value, const std::string attrname);

  // Write node data to HDF5 data file.
  void writeNodeDataReal(const Epetra_Vector &x, const std::string varname);
  void writeNodeDataInt(const Epetra_Vector &x, const std::string varname);

  // Write cell data to HDF5 data file.
  void writeCellDataReal(const Epetra_Vector &x, const std::string varname);
  void writeCellDataInt(const Epetra_Vector &x, const std::string varname);
  
  // Write array data to HDF5 data file. Meant for Restart ONLY not Viz!
  void writeDataReal(const Epetra_Vector &x, const std::string varname);
  void writeDataInt(const Epetra_Vector &x, const std::string varname);
  
  // Read array data from HDF5 data file.
  void readData(Epetra_Vector &x, const std::string varname);
  
 private:

  void createXdmfMesh_(const std::string filename);
  void createXdmfVisit_();

  Teuchos::XMLObject addXdmfHeaderGlobal_();
  Teuchos::XMLObject addXdmfHeaderLocal_(const double value);
  Teuchos::XMLObject addXdmfTopo_();
  Teuchos::XMLObject addXdmfGeo_();
  Teuchos::XMLObject addXdmfAttribute_(std::string varname,
                                       std::string location,
                                       int length, std::string h5path);

  Teuchos::XMLObject findGridNode_(Teuchos::XMLObject xmlobject);
  Teuchos::XMLObject findMeshNode_(Teuchos::XMLObject xmlobject);
  void writeXdmfVisitGrid_(std::string filename);

  void writeFieldData_(const Epetra_Vector &x, std::string varname,
                       datatype_t type, std::string loc);
  void readFieldData_(Epetra_Vector &x, std::string varname,
                      datatype_t type);

  int getCellTypeID_(int conn_len);
  
  
  std::string stripFilename_(std::string filename);
  
  // parallel info
  Epetra_MpiComm viz_comm_;
  MPI_Info info_;
  iogroup_conf_t IOconfig_;
  iogroup_t IOgroup_;
  
  // track xml for viz
  // TODO(barker): need to set default value
  bool TrackXdmf_;

  // XMLObjects for Xdmf output
  Teuchos::XMLObject xmlVisit_;
  Teuchos::XMLObject xmlStep_;

  // Filenames
  std::string H5MeshFilename_;
  std::string H5DataFilename_;
  std::string xdmfVisitFilename_;
  std::string xdmfStepFilename_;

  // Simulation/Mesh Info
  int NumNodes_;
  int NumElems_;
  int ConnLength_;
  int Iteration_;

  // Mesh information
  int conn_;
  AmanziMesh::Cell_type ctype_;
  std::string cname_;
  
  static std::string xdmfHeader_;
};
  
} // close namespace HDF5

#endif  // HDF5MPI_MESH_HH_
