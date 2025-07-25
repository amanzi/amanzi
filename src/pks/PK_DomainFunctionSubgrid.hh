/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

/*
  Process Kernels

*/

#ifndef AMANZI_PK_DOMAIN_FUNCTION_SUBGRID_HH_
#define AMANZI_PK_DOMAIN_FUNCTION_SUBGRID_HH_

#include <string>
#include <vector>

#include "Epetra_Vector.h"
#include "Teuchos_RCP.hpp"

#include "CommonDefs.hh"
#include "DenseVector.hh"
#include "Mesh.hh"
#include "State.hh"
#include "Tag.hh"


namespace Amanzi {

template<class FunctionBase>
class PK_DomainFunctionSubgrid : public FunctionBase {
 public:
  PK_DomainFunctionSubgrid(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh)
    : mesh_(mesh) {};
  virtual ~PK_DomainFunctionSubgrid() = default;

  // member functions
  void Init(const Teuchos::ParameterList& plist,
            const std::string& keyword,
            AmanziMesh::Entity_kind kind);

  // required member functions
  virtual void Compute(double t0, double t1) override;
  virtual DomainFunction_kind getType() const override { return DomainFunction_kind::SUBGRID; }

  virtual void set_state(const Teuchos::RCP<State>& S) final { S_ = S; }

 protected:
  using FunctionBase::keyword_;
  using FunctionBase::value_;

  Teuchos::RCP<const State> S_;
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;

 private:
  AmanziMesh::Entity_ID entity_lid_;
  AmanziMesh::Entity_ID entity_gid_out_;

  Key disaggregated_key_;
  Tag disaggregated_tag_;
};


/* ******************************************************************
* Initialization adds a single function to the list of unique specs.
****************************************************************** */
template<class FunctionBase>
void
PK_DomainFunctionSubgrid<FunctionBase>::Init(const Teuchos::ParameterList& plist,
                                             const std::string& keyword,
                                             AmanziMesh::Entity_kind region_kind)
{
  keyword_ = keyword;

  Teuchos::ParameterList blist = plist.sublist(keyword);
  Key domain_name = Keys::readDomain(blist, "domain", "domain");
  disaggregated_key_ = Keys::readKey(blist, domain_name, "disaggregated");
  disaggregated_tag_ = Keys::readTag(blist, "disaggregated");

  auto regions = plist.get<Teuchos::Array<std::string>>("regions");
  if (regions.size() != 1) {
    Errors::Message m;
    m << "Invalid \"regions\" list in processing subgrid source spec (must be of length 1).";
    Exceptions::amanzi_throw(m);
  }
  std::string region = regions[0];
  if (mesh_->isValidSetName(region, region_kind)) {
    auto id_list = mesh_->getSetEntities(region, region_kind, AmanziMesh::Parallel_kind::ALL);

    if (id_list.size() != 1) {
      Errors::Message m;
      m << "Invalid region in processing subgrid source spec (must be of length 1): \"" << region
        << "\"";
      Exceptions::amanzi_throw(m);
    }
    entity_lid_ = id_list[0];

  } else {
    Errors::Message m;
    m << "Unknown region in processing subgrid source spec: \"" << region << "\"";
    Exceptions::amanzi_throw(m);
  }

  entity_gid_out_ = plist.get<int>("entity GID");
}


/* ******************************************************************
* Compute and distribute the result by Subgrid.
****************************************************************** */
template<class FunctionBase>
void
PK_DomainFunctionSubgrid<FunctionBase>::Compute(double t0, double t1)
{
  auto& vec_out = S_->Get<CompositeVector>(disaggregated_key_, disaggregated_tag_);
  const auto& disaggregated = *vec_out.ViewComponent("cell", true);

  int num_vec = disaggregated.NumVectors();
  std::vector<double> val(num_vec);

  AmanziMesh::Entity_ID entity_lid_out =
    vec_out.Mesh()->getMap(AmanziMesh::Entity_kind::CELL, false).LID(entity_gid_out_);
  AMANZI_ASSERT(entity_lid_out >= 0);
  for (int k = 0; k < num_vec; ++k) val[k] = disaggregated[k][entity_lid_out];
  value_[entity_lid_] = val;
}

} // namespace Amanzi

#endif
