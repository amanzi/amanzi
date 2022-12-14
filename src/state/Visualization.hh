/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Markus Berndt
           Ethan Coon (ecoon@lanl.gov)
*/

//! Manages simulation output to disk.
/*!

A user may request periodic writes of field data for the purposes of
visualization in the `"visualization`" sublists.

ATS accepts a visualization list for each domain/mesh, including surface and
column meshes.  These are in separate ParameterLists, entitled
`"visualization`" for the main mesh, and `"visualization surface`" on the
surface mesh.  It is expected that, for any addition meshes, each will have a
domain name and therefore admit a spec of the form: `"visualization
DOMAIN-NAME`".

.. _visualization-spec:
.. admonition:: visualization-spec

    * `"file name base`" ``[string]`` **visdump_DOMAIN_data**
    * `"dynamic mesh`" ``[bool]`` **false** Write mesh data for every
      visualization dump; this facilitates visualizing deforming meshes.
    * `"time unit`" ``[string]`` **s** A valid time unit to convert time
      into for output files.  One of `"s`", `"d`", `"y`", or `"yr 365`"

    INCLUDES:
    - ``[io-event-spec]`` An IOEvent_ spec


Example:

.. code-block:: xml

  <ParameterList name="visualization">
    <Parameter name="file name base" type="string" value="visdump_data"/>

    <Parameter name="cycles start period stop" type="Array(int)" value="{0, 100, -1}" />
    <Parameter name="cycles" type="Array(int)" value="{999, 1001}" />

    <Parameter name="times start period stop 0" type="Array(double)" value="{0.0, 10.0, 100.0}"/>
    <Parameter name="times start period stop 1" type="Array(double)" value="{100.0, 25.0, -1.0}"/>
    <Parameter name="times" type="Array(double)" value="{101.0, 303.0, 422.0}"/>

    <Parameter name="dynamic mesh" type="bool" value="false"/>
  </ParameterList>

*/

#ifndef AMANZI_STATE_VISUALIZATION_HH_
#define AMANZI_STATE_VISUALIZATION_HH_

#include <vector>
#include <string>

#include "Teuchos_ParameterList.hpp"
#include "Epetra_MultiVector.h"

#include "Units.hh"
#include "Mesh.hh"

#include "IOEvent.hh"
#include "Output.hh"
#include "Tag.hh"

namespace Amanzi {

class Visualization : public IOEvent {
 public:
  Visualization(Teuchos::ParameterList& plist);
  Visualization();

  Teuchos::RCP<const AmanziMesh::Mesh> mesh() const { return mesh_; }
  void set_mesh(const Teuchos::RCP<const AmanziMesh::Mesh> mesh) { mesh_ = mesh; }

  std::string get_name() const { return name_; }
  void set_name(const std::string& name);
  void AddDomain(const std::string& name);
  bool WritesDomain(const std::string& name) const;

  Tag get_tag() const { return tag_; }
  void set_tag(const Tag& tag) { tag_ = tag; }

  // public interface for coordinator clients
  void CreateFiles(bool include_io_set = true);
  void CreateTimestep(double time, int cycle, const std::string& tag);
  virtual void FinalizeTimestep() const;

  // public interface for data clients
  template <typename T>
  void Write(const std::string& name, const T& t) const;

  virtual void
  WriteVector(const Epetra_MultiVector& vec, const std::vector<std::string>& names) const;
  virtual void WriteVector(const Epetra_Vector& vec, const std::string& name) const;
  virtual void WriteRegions();
  virtual void WritePartition();

 protected:
  virtual void ReadParameters_();

  std::vector<std::string> domains_;
  std::string my_units_;
  std::string name_;
  Tag tag_;
  bool time_unit_written_;

  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  Teuchos::RCP<Output> visualization_output_;

  std::map<std::string, Teuchos::Array<std::string>> regions_;
  bool write_partition_;
  bool dynamic_mesh_;
  bool write_mesh_exo_;
};


template <>
inline void
Visualization::Write<Epetra_Vector>(const std::string& name, const Epetra_Vector& t) const
{
  WriteVector(t, name);
}

template <>
inline void
Visualization::Write<double>(const std::string& name, const double& t) const
{
  visualization_output_->WriteAttribute(t, name);
}

template <>
inline void
Visualization::Write<int>(const std::string& name, const int& t) const
{
  visualization_output_->WriteAttribute(t, name);
}

} // namespace Amanzi

#endif
