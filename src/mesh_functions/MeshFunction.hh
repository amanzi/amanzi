/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

//! <MISSING_ONELINE_DOCSTRING>
#pragma once

#include <utility>
#include <vector>
#include <string>

#include "Teuchos_RCP.hpp"

#include "Key.hh"
#include "Mesh.hh"
#include "MultiFunction.hh"
#include "Patch.hh"
#include "CompositeVectorSpace.hh"
#include "CompositeVector.hh"

namespace Amanzi {
namespace Functions {


//
// Class used to hold space and a functor to evaluate on that space
//
class MeshFunction {

public:
  using Spec = std::tuple<std::string, Teuchos::RCP<PatchSpace>, Teuchos::RCP<const MultiFunction>>;
  using SpecList = std::vector<Spec>;

  MeshFunction(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
               AmanziMesh::Entity_kind entity_kind = AmanziMesh::Entity_kind::UNKNOWN);

  MeshFunction(Teuchos::ParameterList& list,
               const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
               const std::string& function_name = "function",
               AmanziMesh::Entity_kind entity_kind = AmanziMesh::Entity_kind::UNKNOWN,
               int flag = 0);

  Teuchos::RCP<const AmanziMesh::Mesh> getMesh() const { return mesh_; }
  void setMesh(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh)
  {
    mesh_ = mesh;
    for (auto& spec : *this) std::get<1>(spec)->mesh = mesh;
  }

  int getFlag() const { return flag_; }

  // add a spec -- others may inherit this and overload to do some checking?
  virtual void addSpec(const Spec& spec);

  void addSpec(const std::string& compname,
               AmanziMesh::Entity_kind entity_kind,
               int num_vectors,
               const std::string& region,
               const Teuchos::RCP<const MultiFunction>& func)
  {
    addSpec(Spec(compname,
                 Teuchos::rcp(new PatchSpace(mesh_, false, region, entity_kind, num_vectors, flag_)),
                 func));
  }

  // access specs
  using spec_iterator = typename SpecList::const_iterator;
  using size_type = typename SpecList::size_type;
  spec_iterator begin() const { return spec_list_.begin(); }
  spec_iterator end() const { return spec_list_.end(); }
  size_type size() const { return spec_list_.size(); }

  using nc_spec_iterator = typename SpecList::iterator;
  nc_spec_iterator begin() { return spec_list_.begin(); }
  nc_spec_iterator end() { return spec_list_.end(); }

  // data creation
  Teuchos::RCP<CompositeVectorSpace> createCVS(bool ghosted) const;
  Teuchos::RCP<MultiPatchSpace> createMPS(bool ghosted = false) const;

  void Compute(double time, MultiPatch<double>& mp);

 protected:
  void readSpec_(Teuchos::ParameterList& list,
                 const std::string& function_name);

 protected:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  SpecList spec_list_;
  AmanziMesh::Entity_kind entity_kind_;
  int flag_;
};


namespace Impl {

//
// helper function to get coordinates, txyz
//
Kokkos::View<double**> getMeshFunctionCoordinates(double time, const PatchSpace& ps);

//
// Computes function on a patch.
//
void
computeFunction(const MultiFunction& f, double time, Patch<double>& p);

//
// Computes function on a patch space, sticking the answer directly into a vector
//
void
computeFunction(const MultiFunction& f, double time, const PatchSpace& p, CompositeVector& cv);

//
// Computes function on a patch space, sticking the answer directly into a vector
//
void
copyFlags(const PatchSpace& p, CompositeVector_<int>& flag_vec);

void
copyFlags(const MultiPatchSpace& mp, CompositeVector_<int>& flag_vec);

void
computeFunctionDepthCoordinate(const MultiFunction& f, double time, Patch<double>& p);


} // namespace Impl


} // namespace Functions
} // namespace Amanzi
