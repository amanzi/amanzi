/*
  State

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Markus Berndt
           Ethan Coon (ecoon@lanl.gov)

  Visualization of data.
*/

#ifndef AMANZI_STATE_VISUALIZATION_HH_
#define AMANZI_STATE_VISUALIZATION_HH_

#include <vector>
#include <string>

#include "Teuchos_ParameterList.hpp"
#include "Epetra_Comm.h"
#include "Epetra_MultiVector.h"

#include "Mesh.hh"

#include "IOEvent.hh"

namespace Amanzi {

class Output;

class Visualization : public IOEvent {
 public:
  Visualization(Teuchos::ParameterList& plist, Epetra_MpiComm *comm);
  Visualization();

  Teuchos::RCP<const AmanziMesh::Mesh> mesh() const { return mesh_; }
  void set_mesh(const Teuchos::RCP<const AmanziMesh::Mesh> mesh) {
    mesh_ = mesh;
  }

  // public interface for coordinator clients
  void CreateFiles();
  void CreateTimestep(const double& time, const int& cycle);
  void FinalizeTimestep() const;

  // public interface for data clients
  void WriteVector(const Epetra_MultiVector& vec, const std::vector<std::string>& names ) const;
  void WriteVector(const Epetra_Vector& vec, const std::string& name ) const;
  void WriteRegions();
  void WritePartition();

 protected:
  void ReadParameters_();

  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  Teuchos::RCP<Output> visualization_output_;

  std::map<std::string, Teuchos::Array<std::string> > regions_;
  bool write_partition_;
  bool dynamic_mesh_;
};

} // Amanzi namespace

#endif
