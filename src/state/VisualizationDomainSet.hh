/*
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Markus Berndt
           Ethan Coon (ecoon@lanl.gov)
*/
//! Visualizes a lifted domain set on the parent mesh.
/*!

Domain sets are collections of subdomains that, most commonly, partition a
paraent domain.  Visualizing this collection of objects can be tricky, but is
made much easier by first migrating all the domain sets on to the parent mesh
and visualizing on that mesh.

.. _visualization-domain-set-spec:
.. admonition:: visualization-domain-set-spec

    INCLUDES:
    - ``[visualization-spec]`` A Visualization_ spec


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

#pragma once

#include <vector>
#include <string>

#include "Teuchos_ParameterList.hpp"

#include "Visualization.hh"

namespace Amanzi {

class VisualizationDomainSet : public Visualization {
 public:
  VisualizationDomainSet(Teuchos::ParameterList& plist) : Visualization(plist) {
    write_partition_ = false; // doesn't work yet
  }

  void set_subdomain_mesh(const Key& subdomain, const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) {
    subdomains_[subdomain] = mesh;
  }

  // public interface for data clients
  virtual void WriteVector(const Epetra_MultiVector& vec, const std::vector<std::string>& names ) const override;
  virtual void WriteVector(const Epetra_Vector& vec, const std::string& name ) const override;

  virtual void FinalizeTimestep() const override;

 protected:
  // note this is lazily constructed, so must be mutable
  mutable std::map<std::string, std::pair<Teuchos::RCP<Epetra_MultiVector>, std::vector<std::string>>> lifted_vectors_;

  std::map<Key, Teuchos::RCP<const AmanziMesh::Mesh>> subdomains_;
};

} // Amanzi namespace

