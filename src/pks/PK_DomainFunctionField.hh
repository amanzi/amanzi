/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatskiy
*/

//! A domain function based on a CompositeVector + Evaluator.
/*
  Process Kernels

*/
/*!

This leverages the Arcos DAG to evaluate a user-provided field and use that in
the domain function.

.. _domain-function-field-spec
.. admonition:: domain-function-field-spec

   * `"field key`" ``[string]`` Field used in the domain function.
   * `"component`" ``[string]`` Component of the field. Default is `"cell`".
   * `"scaling factor`" ``[double]`` Constant multiplication factor. Default is 1.

*/

#ifndef AMANZI_PK_DOMAIN_FUNCTION_FIELD_HH_
#define AMANZI_PK_DOMAIN_FUNCTION_FIELD_HH_

#include <string>
#include <vector>

#include "Epetra_Vector.h"
#include "Teuchos_RCP.hpp"

#include "CommonDefs.hh"
#include "Mesh.hh"
#include "Evaluator.hh"

#include "PK_Utils.hh"
#include "PKsDefs.hh"

namespace Amanzi {

template <class FunctionBase>
class PK_DomainFunctionField : public FunctionBase {
 public:
  PK_DomainFunctionField(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                         const Teuchos::RCP<State>& S,
                         AmanziMesh::Entity_kind kind)
    : S_(S), mesh_(mesh), kind_(kind){};

  PK_DomainFunctionField(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                         const Teuchos::RCP<State>& S,
                         const Teuchos::ParameterList& plist,
                         AmanziMesh::Entity_kind kind)
    : S_(S), mesh_(mesh), FunctionBase(plist), kind_(kind){};

  ~PK_DomainFunctionField(){};

  typedef std::vector<std::string> RegionList;
  typedef std::pair<RegionList, AmanziMesh::Entity_kind> Domain;
  typedef std::set<AmanziMesh::Entity_ID> MeshIDs;


  // member functions
  void Init(const Teuchos::ParameterList& plist, const std::string& keyword);

  // required member functions
  virtual void Compute(double t0, double t1) override;
  virtual DomainFunction_kind getType() const override { return DomainFunction_kind::FIELD; }

 protected:
  using FunctionBase::value_;
  using FunctionBase::keyword_;
  Teuchos::RCP<State> S_;
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;

 private:
  std::string submodel_;
  Teuchos::RCP<MeshIDs> entity_ids_;
  AmanziMesh::Entity_kind kind_;
  Key component_key_;
  int num_cells_;
  std::vector<Key> field_keys_;
  std::vector<Tag> tags_;
};


/* ******************************************************************
* Initialization adds a single function to the list of unique specs.
****************************************************************** */
template <class FunctionBase>
void
PK_DomainFunctionField<FunctionBase>::Init(const Teuchos::ParameterList& plist,
                                           const std::string& keyword)
{
  Errors::Message msg;
  keyword_ = keyword;
  submodel_ = "rate";

  if (plist.isParameter("submodel")) submodel_ = plist.get<std::string>("submodel");
  std::vector<std::string> regions = plist.get<Teuchos::Array<std::string>>("regions").toVector();

  try {
    Teuchos::ParameterList flist = plist.sublist(keyword);
    if (flist.isParameter("number of fields")) {
      // for multiple (>1) field evaluators. This requires "number of fields" sublist
      if (flist.isType<int>("number of fields")) {
        int num_fields = flist.get<int>("number of fields");
        if (num_fields < 2) { // ERROR -- invalid number of fields
          std::cout << "Number of fields must be > 1" << std::endl << std::flush;
          AMANZI_ASSERT(0);
        }

        // Get corresponding field_keys to global memory
        for (int lcv = 1; lcv != (num_fields + 1); ++lcv) {
          std::stringstream sublist_name;
          sublist_name << "field " << lcv << " info";
          Key field_key_ = flist.sublist(sublist_name.str()).get<std::string>("field key");
          Tag tag = Keys::readTag(flist.sublist(sublist_name.str()), "tag");
          field_keys_.push_back(field_key_);
          tags_.push_back(tag);
        }
      }
    } else { // for only ONE field evaluator (without number of fields sublist)
      Key field_key_ = flist.get<std::string>("field key");
      field_keys_.push_back(field_key_);
      Tag tag = Keys::readTag(flist, "tag");
      tags_.push_back(tag);
      // component_key_ = flist.get<std::string>("component", "cell");
    }
    num_cells_ = mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::ALL);
  } catch (Errors::Message& msg) {
    Errors::Message m;
    m << "error in source sublist : " << msg.what();
    Exceptions::amanzi_throw(m);
  }

  // Add this source specification to the domain function.
  Teuchos::RCP<Domain> domain = Teuchos::rcp(new Domain(regions, kind_));

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
* Update the field and stick it in the value map.
****************************************************************** */
template <class FunctionBase>
void
PK_DomainFunctionField<FunctionBase>::Compute(double t0, double t1)
{
  int num_fields = field_keys_.size();      // number of field evaluators
  std::vector<double> val_vec(num_fields);  // vector values
  std::vector<std::vector<double>> field_all(num_fields, std::vector<double>(num_cells_));

  // For multiple field evaluators (fe), loop through each fe to retrieve the values for each cell.
  // Copy the retrieved values to intermediate variables.
  // Finally, update the source values (value_) for each cell.
  for (size_t ife = 0; ife < num_fields; ++ife) { // fe counter
    Key field_key_ = field_keys_[ife];
    Tag tag = tags_[ife];

    if (S_->HasEvaluator(field_key_, tag)) {
      S_->GetEvaluator(field_key_, tag).Update(*S_, field_key_);
    }
    // reference to coressponding fe
    const Epetra_MultiVector& field_vec = *S_->Get<CompositeVector>(field_key_, tag).ViewComponent("cell", false);

    // loop through every cell to update the intermediate variables
    for (auto c : *entity_ids_) {
      field_all[ife][c] = field_vec[0][c];
    }
  }

  // Update source values (value_)
  for (auto c : *entity_ids_) {
    for (int i = 0; i < num_fields; ++i) { val_vec[i] = field_all[i][c]; }
    value_[c] = val_vec;
  }
}

} // namespace Amanzi

#endif
