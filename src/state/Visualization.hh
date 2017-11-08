/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
//! Visualization: a class for controlling simulation output.
/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Markus Berndt
           Ethan Coon (ecoon@lanl.gov)
*/

/*!

Each list contains all parameters as in a IOEvent_ spec, and also:

* `"file name base`" ``[string]`` **"visdump_data"**, **"visdump_surface_data"**
  
* `"dynamic mesh`" ``[bool]`` **false**

  Write mesh data for every visualization dump, this facilitates visualizing deforming meshes.


Example:

.. code-block:: xml

  <ParameterList name="visualization">
    <Parameter name="file name base" type="string" value="visdump_data"/>
  
    <Parameter name="cycles start period stop" type="Array(int)" value="{{0, 100, -1}}" />
    <Parameter name="cycles" type="Array(int)" value="{{999, 1001}}" />

    <Parameter name="times start period stop 0" type="Array(double)" value="{{0.0, 10.0, 100.0}}"/>
    <Parameter name="times start period stop 1" type="Array(double)" value="{{100.0, 25.0, -1.0}}"/>
    <Parameter name="times" type="Array(double)" value="{{101.0, 303.0, 422.0}}"/>

    <Parameter name="dynamic mesh" type="bool" value="false"/>
  </ParameterList>

*/

#ifndef AMANZI_STATE_VISUALIZATION_HH_
#define AMANZI_STATE_VISUALIZATION_HH_

#include <vector>
#include <string>

#include "Teuchos_ParameterList.hpp"
#include "Epetra_MultiVector.h"

#include "Mesh.hh"

#include "IOEvent.hh"

namespace Amanzi {

class Output;

class Visualization : public IOEvent {
 public:
  Visualization(Teuchos::ParameterList& plist);
  Visualization();

  Teuchos::RCP<const AmanziMesh::Mesh> mesh() const { return mesh_; }
  void set_mesh(const Teuchos::RCP<const AmanziMesh::Mesh> mesh) {
    mesh_ = mesh;
  }

  std::string name() const { return name_; }
  void set_name(const std::string& name) { name_ = name; }

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

  std::string name_;
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  Teuchos::RCP<Output> visualization_output_;

  std::map<std::string, Teuchos::Array<std::string> > regions_;
  bool write_partition_;
  bool dynamic_mesh_;
};

} // Amanzi namespace

#endif
