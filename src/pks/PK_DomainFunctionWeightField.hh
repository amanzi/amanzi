/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*!

The total source Q is given for each domain. A weighted source
distribution model is employed. The cell-based source density is
calculated as (Q / W_D) * weight, where W_D is the weighted
domain volume. The weight is defined globally, for the whole
computational domain.

*/

#ifndef AMANZI_PK_DOMAIN_FUNCTION_WEIGHT_FIELD_HH_
#define AMANZI_PK_DOMAIN_FUNCTION_WEIGHT_FIELD_HH_

#include <string>
#include <vector>

#include "Epetra_Vector.h"
#include "Teuchos_RCP.hpp"

#include "Mesh.hh"
#include "Mesh_Algorithms.hh"
#include "State.hh"
#include "Tag.hh"

#include "PK_DomainFunctionWeight.hh"

namespace Amanzi {

template <class FunctionBase>
class PK_DomainFunctionWeightField : public PK_DomainFunctionWeight<FunctionBase> {
 public:
  PK_DomainFunctionWeightField(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                               const Teuchos::RCP<State>& S,
                               AmanziMesh::Entity_kind kind)
    : PK_DomainFunctionWeight<FunctionBase>(mesh, kind), mesh_(mesh), S_(S), kind_(kind) {};
  ~PK_DomainFunctionWeightField(){};

  // member functions
  void Init(const Teuchos::ParameterList& plist, const std::string& keyword);

  // required member functions
  virtual void Compute(double t0, double t1);
  virtual std::string name() const { return "weight by field"; }

 protected:
  using FunctionBase::value_;
  using FunctionBase::domain_volume_;
  using FunctionBase::keyword_;

 private:
  Teuchos::RCP<State> S_;
  Key field_key_;
  Tag tag_;

  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  AmanziMesh::Entity_kind kind_;
};


/* ******************************************************************
* Initialization adds a single function to the list of unique specs.
****************************************************************** */
template <class FunctionBase>
void
PK_DomainFunctionWeightField<FunctionBase>::Init(const Teuchos::ParameterList& plist,
                                                 const std::string& keyword)
{
  field_key_ = plist.get<std::string>("field key");
  tag_ = Tags::DEFAULT;

  std::string name = to_string(kind_); 
  if (S_->Get<CompositeVector>(field_key_, tag_).HasComponent(name)) {
    auto weight = S_->Get<CompositeVector>(field_key_, tag_).ViewComponent(name, true);
    PK_DomainFunctionWeight<FunctionBase>::Init(plist, keyword, weight);
  } else if (kind_ == AmanziMesh::Entity_kind::FACE &&
             S_->Get<CompositeVector>(field_key_, tag_).HasComponent("cell")) {

    auto weight = Teuchos::rcp(new Epetra_MultiVector(mesh_->getMap(AmanziMesh::Entity_kind::FACE, true), 1));
    const auto& field_c = *S_->Get<CompositeVector>(field_key_, tag_).ViewComponent("cell", true);

    int nbf_wghost = mesh_->getNumEntities(AmanziMesh::Entity_kind::BOUNDARY_FACE, AmanziMesh::Parallel_kind::ALL);
    for (int bf = 0; bf < nbf_wghost; ++bf) {
      int f = getBoundaryFaceFace(*mesh_, bf);
      int c = getFaceOnBoundaryInternalCell(*mesh_, f);
      (*weight)[0][f] = field_c[0][c];
    }  
    PK_DomainFunctionWeight<FunctionBase>::Init(plist, keyword, weight);
  }
}


/* ******************************************************************
* Compute and distribute the result by volume.
****************************************************************** */
template <class FunctionBase>
void
PK_DomainFunctionWeightField<FunctionBase>::Compute(double t0, double t1)
{
  // update weight_

  PK_DomainFunctionWeight<FunctionBase>::Compute(t0, t1);
}

} // namespace Amanzi

#endif
