/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*!

This class populates region's data on a derived mesh using a field defined on
the parent mesh. It extends capabilities provided by the "coupling" class by
(a) including interpolation on the parent mesh, (b) rigorous treatment of 
of special cases such as master child entity -> ghost parent entity, which 
may happen for nodes and edges of the derived mesh. 

*/

#ifndef AMANZI_PK_DOMAIN_FUNCTION_PARENT_MESH_FIELD_HH_
#define AMANZI_PK_DOMAIN_FUNCTION_PARENT_MESH_FIELD_HH_

#include <string>
#include <vector>

// TPLs
#include "Epetra_Vector.h"
#include "Teuchos_RCP.hpp"

// Amanzi
#include "CommonDefs.hh"
#include "DenseVector.hh"
#include "Mesh.hh"
#include "MeshHelpers.hh"
#include "PKsDefs.hh"
#include "State.hh"
#include "UniqueLocalIndex.hh"

namespace Amanzi {

template <class FunctionBase>
class PK_DomainFunctionParentMeshField : public FunctionBase {
 public:
  PK_DomainFunctionParentMeshField(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                                   const Teuchos::RCP<const State>& S)
    : mesh_(mesh), S_(S){};
  virtual ~PK_DomainFunctionParentMeshField() = default;

  typedef std::vector<std::string> RegionList;
  typedef std::pair<RegionList, AmanziMesh::Entity_kind> Domain;
  typedef std::set<AmanziMesh::Entity_ID> MeshIDs;

  // member functions
  void Init(const Teuchos::ParameterList& plist,
            const std::string& keyword,
            AmanziMesh::Entity_kind kind);

  // required member functions
  virtual void Compute(double t0, double t1) override;
  virtual DomainFunction_kind getType() const override { return DomainFunction_kind::PARENT_MESH_FIELD; }

 protected:
  using FunctionBase::value_;
  using FunctionBase::keyword_;

  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  Teuchos::RCP<const State> S_;

 private:
  Key field_key_;
  Tag field_tag_;

  Teuchos::RCP<MeshIDs> entity_ids_;
  std::map<AmanziMesh::Entity_ID, AmanziMesh::Entity_ID> reverse_map_;
};


/* ******************************************************************
* Initialization adds a single function to the list of unique specs.
****************************************************************** */
template <class FunctionBase>
void
PK_DomainFunctionParentMeshField<FunctionBase>::Init(const Teuchos::ParameterList& plist,
                                                     const std::string& keyword,
                                                     AmanziMesh::Entity_kind region_kind)
{
  Errors::Message msg;

  keyword_ = keyword;

  Teuchos::ParameterList slist;
  try {
    slist = plist.sublist(keyword_);
  } catch (Errors::Message& msg2) {
    msg << "error in domain coupling sublist : " << msg2.what();
    Exceptions::amanzi_throw(msg);
  }

  // FIXME
  AMANZI_ASSERT(region_kind == AmanziMesh::Entity_kind::FACE);

  // get field name and tag
  field_key_ = slist.get<std::string>("external field key");
  field_tag_ = Keys::readTag(slist, "external field tag");

  // create a list of domain ids
  RegionList regions = plist.get<Teuchos::Array<std::string>>("regions").toVector();
  Teuchos::RCP<Domain> domain = Teuchos::rcp(new Domain(regions, region_kind));

  entity_ids_ = Teuchos::rcp(new MeshIDs());
  AmanziMesh::Entity_kind kind = domain->second;

  for (auto region = domain->first.begin(); region != domain->first.end(); ++region) {
    if (mesh_->isValidSetName(*region, kind)) {
      auto id_list = mesh_->getSetEntities(*region, kind, AmanziMesh::Parallel_kind::OWNED);
      entity_ids_->insert(id_list.begin(), id_list.end());
    } else {
      msg << "Unknown region in processing coupling source: name=" << *region << ", kind=" << kind
          << "\n";
      Exceptions::amanzi_throw(msg);
    }
  }
}


/* ******************************************************************
* Compute and distribute the result by coupling.
****************************************************************** */
template <class FunctionBase>
void
PK_DomainFunctionParentMeshField<FunctionBase>::Compute(double t0, double t1)
{
  auto& field = S_->Get<CompositeVector>(field_key_, field_tag_);
  field.ScatterMasterToGhosted("face");

  auto parent_mesh = field.Mesh();
  const auto& field_f = *field.ViewComponent("face", true);
  auto& fmap = field_f.Map();

  int num_vec = field_f.NumVectors();
  std::vector<double> val(num_vec);

  // Loop over faces (owned + ghosted) on the manifold.
  // Since parent of a master face could be a ghost edge, we work with attached cells. FIXME
  for (auto id = entity_ids_->begin(); id != entity_ids_->end(); ++id) {
    int c = getFaceOnBoundaryInternalCell(*mesh_, *id);
    int f = mesh_->getEntityParent(AmanziMesh::Entity_kind::CELL, c);
    int g = fmap.FirstPointInElement(f);

    for (int k = 0; k < num_vec; ++k) val[k] = (field_f[k][g] + field_f[k][g + 1]) / 2;
    value_[*id] = val;
  }
}

} // namespace Amanzi

#endif
