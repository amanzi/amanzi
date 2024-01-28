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

#include "Units.hh"
#include "Mesh.hh"

#include "IOEvent.hh"
#include "Output.hh"
#include "Tag.hh"
#include "Key.hh"

namespace Amanzi {

class State;

class Visualization : public IOEvent {
 public:
  Visualization(Teuchos::ParameterList& plist,
                const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                bool include_io_set);
  Visualization();

  void addDomain(const std::string& name);
  bool writesDomain(const std::string& name) const;

  // public interface for coordinator clients
  void createFiles(bool include_io_set = true);
  void createTimestep(double time, int cycle);
  virtual void finalizeTimestep();

  // public interface for data clients
  template <typename T>
  void write(const Teuchos::ParameterList& attrs, const T& t) const {
    output_->write(attrs, t);
  }

  // a few special cases get intercepted and delegated to a virtual protected member
  inline void write(const Teuchos::ParameterList& attrs, const Vector_type& t) const {
    writeVector_(attrs, t);
  }
  inline void write(const Teuchos::ParameterList& attrs, const MultiVector_type& t) const {
    writeVector_(attrs, t);
  }

  void writeRegions() const;;
  void writePartition() const;

  void write(const State& S);

 protected:
  void readParameters_();

  virtual void writeVector_(const Teuchos::ParameterList& attrs, const Vector_type& t) const {
    output_->write(attrs, t);
  }
  virtual void writeVector_(const Teuchos::ParameterList& attrs, const MultiVector_type& t) const {
    output_->write(attrs, t);
  }

 protected:
  std::vector<std::string> domains_;
  std::string my_units_;
  std::string name_;
  bool time_unit_written_;
  int count_;

  std::unique_ptr<Output> output_;
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;

  bool write_partition_;
  bool dynamic_mesh_;
  bool write_mesh_exo_;
  std::map<std::string, Teuchos::Array<std::string>> regions_;
};


} // namespace Amanzi

#endif
