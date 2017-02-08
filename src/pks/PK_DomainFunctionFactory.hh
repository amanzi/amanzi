/*
  Process Kernels

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov
*/

#ifndef AMANZI_PK_DOMAIN_FUNCTION_FACTORY_HH_
#define AMANZI_PK_DOMAIN_FUNCTION_FACTORY_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Mesh.hh"
#include "PK_DomainFunctionSimple.hh"
#include "PK_DomainFunctionVolume.hh"
#include "PK_DomainFunctionVolumeFraction.hh"
#include "PK_DomainFunctionWeight.hh"
#include "PK_DomainFunctionCoupling.hh"

namespace Amanzi {


  // Template template class
  template <class ValueType, template <typename Type> class TemplateType>
  class Bar : TemplateType<ValueType>
  {
    TemplateType<ValueType> m_ints;
  };

template <class ValueType, template <typename Type> class FunctionBase>
class PK_DomainFunctionFactory : public FunctionBase<ValueType> {
 public:
  PK_DomainFunctionFactory(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh)
      : mesh_(mesh) {};
  ~PK_DomainFunctionFactory() {};

  Teuchos::RCP<FunctionBase<ValueType> > Create(const Teuchos::ParameterList& plist,
                                                AmanziMesh::Entity_kind kind,
                                                Teuchos::RCP<const Epetra_Vector> weight);

  Teuchos::RCP<FunctionBase<ValueType> > Create(const Teuchos::ParameterList& plist,
                                                const std::string& keyword,
                                                AmanziMesh::Entity_kind kind,
                                                Teuchos::RCP<const Epetra_Vector> weight);

 protected:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
};


template <class ValueType, template <typename Type> class FunctionBase>
Teuchos::RCP<FunctionBase<ValueType> > 
PK_DomainFunctionFactory<ValueType, FunctionBase >::Create(
    const Teuchos::ParameterList& plist,
    const std::string& keyword,
    AmanziMesh::Entity_kind kind,
    Teuchos::RCP<const Epetra_Vector> weight)
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
    Teuchos::RCP<PK_DomainFunctionVolumeFraction<ValueType, FunctionBase> >
      func = Teuchos::rcp(new PK_DomainFunctionVolumeFraction<ValueType, FunctionBase>(mesh_, kind));
    func->Init(plist, keyword);
    return func; 
  }
  else if (model == "volume") {
    Teuchos::RCP<PK_DomainFunctionVolume<ValueType, FunctionBase> >
        func = Teuchos::rcp(new PK_DomainFunctionVolume<ValueType, FunctionBase>(mesh_, plist, kind));
    func->Init(plist, keyword);
    return func; 
  }
  else if (model == "permeability") {
    Teuchos::RCP<PK_DomainFunctionWeight<ValueType, FunctionBase> >
        func = Teuchos::rcp(new PK_DomainFunctionWeight<ValueType, FunctionBase>(mesh_));
    func->Init(plist, keyword, weight);
    return func; 
  }
  else if (model == "domain coupling") {
    Teuchos::RCP<PK_DomainFunctionCoupling<ValueType, FunctionBase> >
      func = Teuchos::rcp(new PK_DomainFunctionCoupling<ValueType, FunctionBase>(mesh_));
    func->Init(plist, keyword, kind);
    return func;
  }
  else {
    Teuchos::RCP<PK_DomainFunctionSimple<ValueType, FunctionBase> >
        func = Teuchos::rcp(new PK_DomainFunctionSimple<ValueType, FunctionBase>(mesh_, plist, kind));
    func->Init(plist, keyword);
    return func; 
  }

  return Teuchos::null;
}


template <class ValueType, template <typename Type>class FunctionBase>
Teuchos::RCP<FunctionBase<ValueType> > PK_DomainFunctionFactory<ValueType, FunctionBase >::Create(
                                               const Teuchos::ParameterList& plist,
                                               AmanziMesh::Entity_kind kind,
                                               Teuchos::RCP<const Epetra_Vector> weight)
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

  return Create(plist, keyword, kind, weight);
}

}  // namespace Amanzi

#endif

