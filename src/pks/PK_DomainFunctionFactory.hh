/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov
*/

/*
  Process Kernels

*/

#ifndef AMANZI_PK_DOMAIN_FUNCTION_FACTORY_HH_
#define AMANZI_PK_DOMAIN_FUNCTION_FACTORY_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Mesh.hh"
#include "State.hh"
#include "PK_DomainFunctionCoupling.hh"
#include "PK_DomainFunctionField.hh"
#include "PK_DomainFunctionFirstOrderExchange.hh"
#include "PK_DomainFunctionSimple.hh"
#include "PK_DomainFunctionSimpleWell.hh"
#include "PK_DomainFunctionSubgrid.hh"
#include "PK_DomainFunctionSubgridReturn.hh"
#include "PK_DomainFunctionVolume.hh"
#include "PK_DomainFunctionVolumeFraction.hh"
#include "PK_DomainFunctionWeightField.hh"
#include "PK_DomainFunctionWeight.hh"

namespace Amanzi {

template <class FunctionBase>
class PK_DomainFunctionFactory : public FunctionBase {
 public:
  PK_DomainFunctionFactory(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                           const Teuchos::RCP<State>& S)
    : mesh_(mesh), S_(S){};
  ~PK_DomainFunctionFactory(){};

  Teuchos::RCP<FunctionBase> Create(Teuchos::ParameterList& plist,
                                    AmanziMesh::Entity_kind kind,
                                    Teuchos::RCP<const Epetra_MultiVector> weight,
                                    const Tag& tag = Tags::DEFAULT);

  Teuchos::RCP<FunctionBase> Create(Teuchos::ParameterList& plist,
                                    const std::string& keyword,
                                    AmanziMesh::Entity_kind kind,
                                    Teuchos::RCP<const Epetra_MultiVector> weight,
                                    const Tag& tag = Tags::DEFAULT);


 protected:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  Teuchos::RCP<State> S_;
};


template <class FunctionBase>
Teuchos::RCP<FunctionBase>
PK_DomainFunctionFactory<FunctionBase>::Create(Teuchos::ParameterList& plist,
                                               const std::string& keyword,
                                               AmanziMesh::Entity_kind kind,
                                               Teuchos::RCP<const Epetra_MultiVector> weight,
                                               const Tag& tag)
{
  // verify completeness of the list
  Errors::Message msg;

  if (!plist.isParameter("regions")) {
    msg << "required parameter \"regions\" is missing";
    Exceptions::amanzi_throw(msg);
  }

  if (!plist.isSublist(keyword)) {
    msg << "required parameter \"" << keyword << "\" is not a sublist";
    Exceptions::amanzi_throw(msg);
  }

  // should we use available volume fractions
  bool use_vofs(false);
  if (plist.isParameter("use area fractions")) {
    use_vofs = plist.get<bool>("use area fractions");
  } else if (plist.isParameter("use volume fractions")) {
    use_vofs = plist.get<bool>("use volume fractions");
  }

  // select model for data distribution
  std::string model = plist.get<std::string>("spatial distribution method");

  if (use_vofs) {
    Teuchos::RCP<PK_DomainFunctionVolumeFraction<FunctionBase>> func =
      Teuchos::rcp(new PK_DomainFunctionVolumeFraction<FunctionBase>(mesh_, kind));
    func->Init(plist, keyword);
    return func;
  } else if (model == "volume") {
    Teuchos::RCP<PK_DomainFunctionVolume<FunctionBase>> func =
      Teuchos::rcp(new PK_DomainFunctionVolume<FunctionBase>(mesh_, plist, kind));
    func->Init(plist, keyword);
    return func;
  } else if (model == "permeability") {
    Teuchos::RCP<PK_DomainFunctionWeight<FunctionBase>> func =
      Teuchos::rcp(new PK_DomainFunctionWeight<FunctionBase>(mesh_, kind));
    func->Init(plist, keyword, weight);
    return func;
  } else if (model == "weight by field") {
    Teuchos::RCP<PK_DomainFunctionWeightField<FunctionBase>> func =
      Teuchos::rcp(new PK_DomainFunctionWeightField<FunctionBase>(mesh_, S_, kind));
    func->Init(plist, keyword);
    return func;
  } else if (model == "domain coupling") {
    plist.set<std::string>("external field copy key", tag.get());
    Teuchos::RCP<PK_DomainFunctionCoupling<FunctionBase>> func =
      Teuchos::rcp(new PK_DomainFunctionCoupling<FunctionBase>(mesh_, S_));
    func->Init(plist, keyword, kind);
    return func;
  } else if (model == "field") {
    Teuchos::RCP<PK_DomainFunctionField<FunctionBase>> func =
      Teuchos::rcp(new PK_DomainFunctionField<FunctionBase>(mesh_, S_, kind));
    func->Init(plist, keyword);
    return func;
  } else if (model == "first order exchange") {
    AMANZI_ASSERT(kind == AmanziMesh::Entity_kind::CELL);
    plist.sublist("source function")
      .set<std::string>("total component concentration copy", tag.get());
    Teuchos::RCP<PK_DomainFunctionFirstOrderExchange<FunctionBase>> func =
      Teuchos::rcp(new PK_DomainFunctionFirstOrderExchange<FunctionBase>(mesh_, plist, kind));
    func->Init(plist, keyword);
    return func;
  } else if (model == "subgrid") {
    AMANZI_ASSERT(kind == AmanziMesh::Entity_kind::FACE);
    plist.sublist(keyword).set<std::string>("copy_field_out_tag", tag.get());
    Teuchos::RCP<PK_DomainFunctionSubgrid<FunctionBase>> func =
      Teuchos::rcp(new PK_DomainFunctionSubgrid<FunctionBase>(mesh_));
    func->Init(plist, keyword, kind);
    return func;
  } else if (model == "subgrid return") {
    AMANZI_ASSERT(kind == AmanziMesh::Entity_kind::CELL);
    plist.sublist("source function").set<std::string>("copy subgrid field", tag.get());
    Teuchos::RCP<PK_DomainFunctionSubgridReturn<FunctionBase>> func =
      Teuchos::rcp(new PK_DomainFunctionSubgridReturn<FunctionBase>(mesh_, plist));
    func->Init(plist, keyword);
    return func;
  } else if (model == "simple well") {
    Teuchos::RCP<PK_DomainFunctionSimpleWell<FunctionBase>> func =
      Teuchos::rcp(new PK_DomainFunctionSimpleWell<FunctionBase>(mesh_));
    func->Init(plist, keyword, S_);
    return func;
  } else {
    Teuchos::RCP<PK_DomainFunctionSimple<FunctionBase>> func =
      Teuchos::rcp(new PK_DomainFunctionSimple<FunctionBase>(mesh_, plist, kind));
    func->Init(plist, keyword);
    return func;
  }

  return Teuchos::null;
}


template <class FunctionBase>
Teuchos::RCP<FunctionBase>
PK_DomainFunctionFactory<FunctionBase>::Create(Teuchos::ParameterList& plist,
                                               AmanziMesh::Entity_kind kind,
                                               Teuchos::RCP<const Epetra_MultiVector> weight,
                                               const Tag& tag)
{
  int n(0);
  std::string keyword;
  for (auto it = plist.begin(); it != plist.end(); ++it) {
    if (plist.isSublist(it->first)) {
      n++;
      keyword = it->first;
    }
  }
  if (keyword.size() == 0 || n > 1) {
    Errors::Message msg;
    msg << "Domain function should have exactly one sublist";
    Exceptions::amanzi_throw(msg);
  }

  return Create(plist, keyword, kind, weight, tag);
}

} // namespace Amanzi

#endif
