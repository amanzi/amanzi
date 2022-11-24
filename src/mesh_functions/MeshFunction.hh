/*
  Mesh Functions

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon

  Function applied to a mesh component.
*/

#ifndef AMANZI_MESH_FUNCTION_HH_
#define AMANZI_MESH_FUNCTION_HH_


#include <utility>
#include <vector>
#include <string>

#include "Teuchos_RCP.hpp"
#include "Mesh.hh"
#include "MultiFunction.hh"

namespace Amanzi {
namespace Functions {

class MeshFunction {
 public:
  // Tuple of a region plus a mesh entity provides the domain on which a
  // function can be evaluated.
  typedef std::vector<std::string> RegionList;
  typedef std::pair<RegionList, AmanziMesh::Entity_kind> Domain;

  // A specification for domain and function.
  typedef std::pair<Teuchos::RCP<Domain>, Teuchos::RCP<const MultiFunction>> Spec;
  typedef std::vector<Teuchos::RCP<Spec>> SpecList;

  // constructor
  MeshFunction(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) : mesh_(mesh){};
  virtual ~MeshFunction() = default;

  // add a spec -- others may inherit this and overload to do some checking?
  virtual void AddSpec(const Teuchos::RCP<Spec>& spec) { spec_list_.push_back(spec); }

  // access specs
  typedef SpecList::const_iterator spec_iterator;
  spec_iterator begin() const { return spec_list_.begin(); }
  spec_iterator end() const { return spec_list_.end(); }
  SpecList::size_type size() const { return spec_list_.size(); }

  // access mesh
  Teuchos::RCP<const AmanziMesh::Mesh> mesh() const { return mesh_; }

 protected:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  SpecList spec_list_;
};

} // namespace Functions
} // namespace Amanzi


#endif
