/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
Amanzi

License:
Author: Markus Berndt
        Ethan Coon (ecoon@lanl.gov)

Visualization of data.

------------------------------------------------------------------------- */

#ifndef AMANZI_STATE_VISUALIZATION_HH_
#define AMANZI_STATE_VISUALIZATION_HH_

#include "Teuchos_ParameterList.hpp"
#include "Epetra_Comm.h"
#include "Mesh.hh"
#include "hdf5mpi_mesh.hh"

#include "io_event.hh"

namespace Amanzi {

class Visualization : public IOEvent {

 public:

  Visualization(Teuchos::ParameterList& plist, Epetra_MpiComm *comm);
  Visualization();

  Teuchos::RCP<const AmanziMesh::Mesh> mesh() const { return mesh_; }
  void set_mesh(const Teuchos::RCP<const AmanziMesh::Mesh> mesh) { mesh_ = mesh; }

  // public interface for coordinator clients
  void CreateFiles();
  void CreateTimestep(const double& time, const int& cycle);
  void FinalizeTimestep() const;

  // public interface for data clients
  void WriteMesh(const double time, const int iteration) const;
  void WriteVector(const Epetra_MultiVector& vec, const std::vector<std::string>& names ) const;
  void WriteVector(const Epetra_Vector& vec, const std::string name ) const;


 protected:
  void ReadParameters_();

  std::string filebasename_;
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  Teuchos::RCP<Amanzi::HDF5_MPI> visualization_output_;

  bool dynamic_mesh_;
}; // Visualization class

} // Amanzi namespace

#endif
