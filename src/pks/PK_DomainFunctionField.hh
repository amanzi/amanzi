/*
  Process Kernels 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_PK_DOMAIN_FUNCTION_FIELD_HH_
#define AMANZI_PK_DOMAIN_FUNCTION_FIELD_HH_

#include <string>
#include <vector>

#include "Epetra_Vector.h"
#include "Teuchos_RCP.hpp"

#include "CommonDefs.hh"
#include "DenseVector.hh"
#include "Mesh.hh"
#include "UniqueMeshFunction.hh"

#include "PK_Utils.hh"

namespace Amanzi {

template <class FunctionBase>
class PK_DomainFunctionField : public FunctionBase,
                                public Functions::UniqueMeshFunction {
 public:
  PK_DomainFunctionField(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                         const Teuchos::RCP<const State>& S,
                          AmanziMesh::Entity_kind kind) :
    S_(S),
    UniqueMeshFunction(mesh),
    kind_(kind) {};

  PK_DomainFunctionField(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                         const Teuchos::RCP<const State>& S,
                         const Teuchos::ParameterList& plist,
                         AmanziMesh::Entity_kind kind) :
    S_(S),
    UniqueMeshFunction(mesh),
    FunctionBase(plist),
    kind_(kind) {
  };

  ~PK_DomainFunctionField() {};

  // member functions
  void Init(const Teuchos::ParameterList& plist, const std::string& keyword);

  // required member functions
  virtual void Compute(double t0, double t1);
  virtual std::string name() const { return "simple"; }

 protected:
  using FunctionBase::value_;
  using FunctionBase::keyword_;  
  Teuchos::RCP<const State> S_;

 private:
  std::string submodel_;
  AmanziMesh::Entity_kind kind_;
  std::string field_key_;
  std::string copy_key_;
  std::string component_key_;
};


/* ******************************************************************
* Initialization adds a single function to the list of unique specs.
****************************************************************** */
template <class FunctionBase>
void PK_DomainFunctionField<FunctionBase>::Init(
    const Teuchos::ParameterList& plist, const std::string& keyword)
{
  keyword_ = keyword;

  submodel_ = "rate";
  if (plist.isParameter("submodel"))
    submodel_ = plist.get<std::string>("submodel");
  std::vector<std::string> regions = plist.get<Teuchos::Array<std::string> >("regions").toVector();
  
  // Teuchos::RCP<Amanzi::MultiFunction> f;
  try {
    Teuchos::ParameterList flist = plist.sublist(keyword);
    field_key_ = flist.get<std::string>("field key");
    copy_key_ = flist.get<std::string>("copy key", "default");
    component_key_ = flist.get<std::string>("component", "cell");;    
  } catch (Errors::Message& msg) {
    Errors::Message m;
    m << "error in source sublist : " << msg.what();
    Exceptions::amanzi_throw(m);
  }

  // Add this source specification to the domain function.
  Teuchos::RCP<Domain> domain = Teuchos::rcp(new Domain(regions, kind_));
  //AddSpec(Teuchos::rcp(new Spec(domain, f)));
}


/* ******************************************************************
* Compute and distribute the result by Simple.
****************************************************************** */
template <class FunctionBase>
void PK_DomainFunctionField<FunctionBase>::Compute(double t0, double t1)
{
  if (unique_specs_.size() == 0) return;

  // create the input tuple (time + space)
  int dim = mesh_->space_dimension();
  const auto& field_vec = *S_->GetFieldCopyData(field_key_, copy_key_)->ViewComponent(component_key_, false);
  int nvalues = field_vec.NumVectors();
  
  for (auto uspec = unique_specs_.at(kind_)->begin(); uspec != unique_specs_.at(kind_)->end(); ++uspec) {

    Teuchos::RCP<MeshIDs> ids = (*uspec)->second;
    // uspec->first is a RCP<Spec>, Spec's second is an RCP to the function.
    //int nfun = (*uspec)->first->second->size();
    std::vector<double> val_vec(nvalues);  

    for (auto c = ids->begin(); c != ids->end(); ++c) {
      for (int i = 0; i < nvalues; ++i) {
        //val_vec[i] = field[i][c];
      }
      value_[*c] = val_vec;
    }

  }
}

}  // namespace Amanzi

#endif
