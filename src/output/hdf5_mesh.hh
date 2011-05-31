#ifndef HDF5_MESH_HH_
#define HDF5_MESH_HH_

#include <string>

#include "Mesh_maps_base.hh"
#include "errors.hh"
#include "Epetra_Vector.h"
#include "Teuchos_XMLObject.hpp"

extern "C" {
#include "hdf5.h"
};

class HDF5 {


 public:

  bool TrackXdmf() { return TrackXdmf_; }
  void setTrackXdmf(bool TrackXdmf) { TrackXdmf_ = TrackXdmf; }

  std::string H5MeshFilename() { return H5MeshFilename_; }
  void setH5MeshFilename(std::string H5MeshFilename) {
                         H5MeshFilename_ = H5MeshFilename;}
  std::string H5DataFilename() { return H5DataFilename_; }
  void setH5DataFilename(std::string H5DataFilename) {
                         H5DataFilename_ = H5DataFilename;}
  std::string xdmfParaviewFilename() { return xdmfParaviewFilename_; }
  void setxdmfParaviewFilename(std::string xdmfParaviewFilename) {
                               xdmfParaviewFilename_ = xdmfParaviewFilename;}
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
  int Iteration() { return Iteration_;}
  void setIteration(int Iteration) {Iteration_ = Iteration;}

  Teuchos::XMLObject xmlParaview() { return xmlParaview_; }
  Teuchos::XMLObject xmlVisit() { return xmlVisit_; }
  Teuchos::XMLObject xmlStep() { return xmlStep_; }

  // Output mesh data to filename.h5 and filename.xmf
  void createMeshFile(Mesh_maps_base &mesh_Maps, std::string filename);

  // Create h5 file for data output, create accompanying Xdmf files for
  // ParaView and Visit
  void createDataFile(std::string mesh_filename, std::string data_filename);
  // Adds time step attributes to ParaView and VisIt Xdmf files.  Creates
  // individual Xdmf for the current step.
  // TODO(barker): The individual step file can be remove after VisIt updates.
  // TODO(barker): Consolidate into a singel Xdmf file, after VisIt updates.
  void createTimestep(const double time, const int iteration);
  void endTimestep();

  // Write node data to HDF5 data file.
  void writeNodeData(const Epetra_Vector &x, const std::string varname);

  // Write cell data to HDF5 data file.
  void writeCellData(const Epetra_Vector &x, const std::string varname);

 private:

  void createXdmfMesh_(const std::string filename);
  void createXdmfParaview_();
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
  void writeXdmfParaviewGrid_(std::string filename,
                                          const double time,
                                          const int iteration);
  void writeXdmfVisitGrid_(std::string filename);

  void writeFieldData_(const Epetra_Vector &x, std::string varname,
                       std::string loc);

  // track xml for viz
  // TODO(barker): need to set default value
  bool TrackXdmf_;

  // XMLObjects for Xdmf output
  Teuchos::XMLObject xmlParaview_;
  Teuchos::XMLObject xmlVisit_;
  Teuchos::XMLObject xmlStep_;

  // Filenames
  std::string H5MeshFilename_;
  std::string H5DataFilename_;
  std::string xdmfParaviewFilename_;
  std::string xdmfVisitFilename_;
  std::string xdmfStepFilename_;

  // Simulation/Mesh Info
  int NumNodes_;
  int NumElems_;
  int Iteration_;

  static std::string xdmfHeader_;
};

#endif  // HDF5_MESH_HH_
