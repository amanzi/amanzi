/*
  Output

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Erin Baker
*/

#ifndef HDF5MPI_MESH_HH_
#define HDF5MPI_MESH_HH_

#include <set>
#include <string>
#include <fstream>

//#ifdef HAVE_MOAB_MESH
//#include "Mesh_moab.hh"
//#endif

#include "Mesh.hh"
#include "MeshDefs.hh"
#include "errors.hh"
#include "dbc.hh"
#include "Epetra_Vector.h"
#include "Teuchos_XMLObject.hpp"

extern "C" {
#include "hdf5.h"
#include "hdf5_hl.h"
#include "parallelIO.h"
};

#define MAX_STRING_LENGTH 100

namespace Amanzi {

class HDF5_MPI {
 public:
  HDF5_MPI(const Comm_ptr_type &comm);
  HDF5_MPI(const Comm_ptr_type &comm, std::string dataFilename);
  ~HDF5_MPI(void);
  
  bool TrackXdmf() { return TrackXdmf_; }
  void setTrackXdmf(bool TrackXdmf) { TrackXdmf_ = TrackXdmf; }

  std::string H5MeshFilename() { return H5MeshFilename_; }
  void setH5MeshFilename(std::string H5MeshFilename) { H5MeshFilename_ = H5MeshFilename;}
  std::string H5DataFilename() { return H5DataFilename_; }
  void setH5DataFilename(std::string H5DataFilename) { H5DataFilename_ = H5DataFilename;}
  std::string xdmfMeshVisitFilename() { return xdmfMeshVisitFilename_; }
  void setxdmfMeshVisitFilename(std::string xdmfMeshVisitFilename) { xdmfMeshVisitFilename_ = xdmfMeshVisitFilename; }  
  std::string xdmfStepFilename() { return xdmfStepFilename_; }
  void setxdmfStepFilename(std::string xdmfStepFilename) { xdmfStepFilename_ = xdmfStepFilename; }
  
  int NumNodes() { return NumNodes_; }
  void setNumNodes(int NumNodes) { NumNodes_ = NumNodes; }
  int NumElems() { return NumElems_; }
  void setNumElems(int NumElems) { NumElems_ = NumElems; }
  int ConnLength() { return ConnLength_; }
  void setConnLength(int ConnLength) { ConnLength_ = ConnLength; }
  int Iteration() { return Iteration_;}
  void setIteration(int Iteration) { Iteration_ = Iteration; }
  double Time() { return Time_;}
  void setTime(double Time) { Time_ = Time; }
  void setDynMesh(const bool dynamic_mesh) { dynamic_mesh_ = dynamic_mesh; }

  Teuchos::XMLObject xmlMeshVisit() { return xmlMeshVisit_; }

  // Output mesh data to filename.h5 and filename.xmf
  void createMeshFile(Teuchos::RCP<const AmanziMesh::Mesh> mesh, const std::string& filename);
  void writeMesh(const double time, const int iteration);
  void writeDualMesh(const double time, const int iteration);

  // Create h5 file for data output, create accompanying Xdmf files for Visit
  void createDataFile(const std::string& data_filename);
  // Adds time step attributes to VisIt Xdmf files.  Creates
  // individual Xdmf for the current step.
  // TODO(barker): Consolidate into a singel Xdmf file, after VisIt updates.
  void createTimestep(const double time, const int iteration);
  void endTimestep();
  
  // before writing data to the h5 file, the user must open the
  // file, and after writing is done, he must close it
  void open_h5file();
  void close_h5file();


  // Write attribute to HDF5 data file.
  void writeAttrReal(double value, const std::string attrname);
  void writeAttrReal(double value, const std::string attrname, std::string h5path);
  void writeAttrReal(double *value, int ndim, const std::string attrname);
  void writeAttrInt(int value, const std::string attrname);
  void writeAttrInt(int *value, int ndim, const std::string attrname);
  void writeAttrString(const std::string value, const std::string attrname);
  void readAttrReal(double &value, const std::string attrname);
  void readAttrReal(double **value, int *ndim, const std::string attrname);
  void readAttrInt(int &value, const std::string attrname);
  void readAttrInt(int **value, int *ndim, const std::string attrname);
  void readAttrString(std::string &value, const std::string attrname);

  // Write node data to HDF5 data file.
  void writeNodeDataReal(const Epetra_Vector &x, const std::string& varname);
  void writeNodeDataInt(const Epetra_Vector &x, const std::string& varname);

  // Write cell data to HDF5 data file.
  void writeCellDataReal(const Epetra_Vector &x, const std::string& varname);
  void writeCellDataInt(const Epetra_Vector &x, const std::string& varname);
  
  // Write array data to HDF5 data file. Meant for Restart ONLY not Viz!
  void writeDataReal(const Epetra_Vector &x, const std::string& varname);
  void writeDataInt(const Epetra_Vector &x, const std::string& varname);
  
  // Read array data from HDF5 data file.
  bool readData(Epetra_Vector &x, const std::string varname);
  
  // Write and read string datasets
  void writeDataString(char **x, int num_entries, const std::string& varname);
  void readDataString(char ***x, int *num_entries, const std::string& varname);
  
  // -- due general lack of parallel distribution, root reads all data
  void writeDatasetReal(double *data, int nloc, int nglb, const std::string& varname);
  bool readDatasetReal(double **data, int nloc, const std::string& varname);

 private:
  void createXdmfMesh_(const std::string filename, const double time, const int iteration);
  void createXdmfVisit_();
  void createXdmfMeshVisit_();

  Teuchos::XMLObject addXdmfHeaderGlobal_();
  Teuchos::XMLObject addXdmfHeaderLocal_(const std::string& name, const double value, const int cycle);
  Teuchos::XMLObject addXdmfTopo_(const int cycle);
  Teuchos::XMLObject addXdmfGeo_(const int cycle);
  Teuchos::XMLObject addXdmfAttribute_(std::string varname,
                                       std::string location,
                                       int length, std::string h5path);

  Teuchos::XMLObject findGridNode_(Teuchos::XMLObject xmlobject);
  Teuchos::XMLObject findMeshNode_(Teuchos::XMLObject xmlobject);
  void writeXdmfVisitGrid_(std::string filename);
  void writeXdmfMeshVisitGrid_(std::string filename);

  void writeFieldData_(const Epetra_Vector &x, const std::string& varname, datatype_t type, std::string loc);
  bool readFieldData_(Epetra_Vector &x, const std::string& varname, datatype_t type);

  bool checkFieldData_(const std::string& varname);

  int getCellTypeID_(AmanziMesh::Cell_type type);
  
  std::set<std::string> extractFields_(const Teuchos::XMLObject& xml);
  
  std::string stripFilename_(std::string filename);
  
  // mesh
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;

  // parallel info
  Teuchos::RCP<const MpiComm_type> viz_comm_;
  MPI_Info info_;
  iogroup_conf_t IOconfig_;
  iogroup_t IOgroup_;
  
 private:
  // track xml for viz
  // TODO(barker): need to set default value
  bool TrackXdmf_;

  // XMLObjects for Xdmf output
  std::vector<Teuchos::XMLObject> xmlVisit_;
  Teuchos::XMLObject xmlMeshVisit_;
  Teuchos::XMLObject xmlStep_, xmlStep_prev_;

  // Filenames
  std::string H5MeshFilename_;
  std::string H5DataFilename_;
  std::string xdmfVisitFilename_, xdmfVisitBaseFilename_;
  std::string xdmfMeshVisitFilename_;
  std::string xdmfStepFilename_;
  std::string base_filename_;
  std::string h5Filename_;

  // Simulation/Mesh Info
  int NumNodes_;
  int NumElems_;
  int ConnLength_;
  int Iteration_;
  double Time_;

  // Mesh information
  std::string cname_;
  
  static std::string xdmfHeader_;

  hid_t mesh_file_;
  hid_t data_file_;
  std::ofstream of_timestep_;

  bool dynamic_mesh_, mesh_written_;
  int static_mesh_cycle_;
};
  
} // close namespace HDF5

#endif  // HDF5MPI_MESH_HH_
