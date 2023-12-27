/*
  Copyright 2010-202x held jointly by participating institutions.
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

    <Parameter name="cycles start period stop" type="Array(int)" value="{0, 100, -1}" />
    <Parameter name="cycles" type="Array(int)" value="{999, 1001}" />

    <Parameter name="times start period stop 0" type="Array(double)" value="{0.0, 10.0, 100.0}"/>
    <Parameter name="times start period stop 1" type="Array(double)" value="{100.0, 25.0, -1.0}"/>
    <Parameter name="times" type="Array(double)" value="{101.0, 303.0, 422.0}"/>

    <Parameter name="dynamic mesh" type="bool" value="false"/>
  </ParameterList>

*/

#pragma once

#include <vector>
#include <string>

#include "Teuchos_ParameterList.hpp"

#include "Visualization.hh"

namespace Amanzi {

namespace AmanziMesh {
class DomainSet;
}

class VisualizationDomainSet : public Visualization {
 public:
  VisualizationDomainSet(Teuchos::ParameterList& plist,
                         const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                         bool include_io_set,
                         const Teuchos::RCP<const AmanziMesh::DomainSet>& ds)
    : Visualization(plist, mesh, include_io_set)
  {
    write_partition_ = false; // doesn't work yet
    setDomainSet(ds);
  }

  void setDomainSet(const Teuchos::RCP<const AmanziMesh::DomainSet>& ds) { ds_ = ds; }

  virtual void finalizeTimestep() override;

 protected:
  // public interface for data clients
  virtual void
  writeVector_(const Teuchos::ParameterList& attrs, const MultiVector_type& vec) const override;
  virtual void
  writeVector_(const Teuchos::ParameterList& attrs, const Vector_type& vec) const override;

 protected:
  // note this is lazily constructed, so must be mutable
  // Note that this relies on map being ORDERED!
  mutable std::map<std::string, std::pair<Teuchos::RCP<MultiVector_type>, Teuchos::ParameterList>>
    lifted_vectors_;

  // this is constructed after all are set, and provides a common set of shared
  // names that exist on any rank for collective writes.
  mutable std::vector<std::string> lifted_vector_names_;
  Teuchos::RCP<const AmanziMesh::DomainSet> ds_;
  std::string dset_name_;
};

} // namespace Amanzi
