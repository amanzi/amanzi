/*
  Process Kernels 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_PK_DOMAIN_FUNCTION_SIMPLE_HH_
#define AMANZI_PK_DOMAIN_FUNCTION_SIMPLE_HH_

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
class PK_DomainFunctionSimple : public FunctionBase,
                                public Functions::UniqueMeshFunction {
 public:
  PK_DomainFunctionSimple(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                          AmanziMesh::Entity_kind kind) :
      UniqueMeshFunction(mesh),
      kind_(kind) {};

  PK_DomainFunctionSimple(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                          const Teuchos::ParameterList& plist,
                          AmanziMesh::Entity_kind kind) :
    FunctionBase(plist),
    UniqueMeshFunction(mesh),
    kind_(kind) {
  };

  ~PK_DomainFunctionSimple() {};

  // member functions
  void Init(const Teuchos::ParameterList& plist, const std::string& keyword);

  // required member functions
  virtual void Compute(double t0, double t1);
  virtual std::string name() const { return "simple"; }

 protected:
  using FunctionBase::value_;
  using FunctionBase::keyword_;

 private:
  std::string submodel_;
  AmanziMesh::Entity_kind kind_;
};


/* ******************************************************************
* Initialization adds a single function to the list of unique specs.
****************************************************************** */
template <class FunctionBase>
void PK_DomainFunctionSimple<FunctionBase>::Init(
    const Teuchos::ParameterList& plist, const std::string& keyword)
{
  keyword_ = keyword;

  submodel_ = "rate";
  if (plist.isParameter("submodel"))
    submodel_ = plist.get<std::string>("submodel");
  std::vector<std::string> regions = plist.get<Teuchos::Array<std::string> >("regions").toVector();

  Teuchos::RCP<Amanzi::MultiFunction> f;
  try {
    Teuchos::ParameterList flist = plist.sublist(keyword);
    f = Teuchos::rcp(new MultiFunction(flist));
  } catch (Errors::Message& msg) {
    Errors::Message m;
    m << "error in source sublist : " << msg.what();
    Exceptions::amanzi_throw(m);
  }

  // Add this source specification to the domain function.
  Teuchos::RCP<Domain> domain = Teuchos::rcp(new Domain(regions, kind_));
  AddSpec(Teuchos::rcp(new Spec(domain, f)));
}


/* ******************************************************************
* Compute and distribute the result by Simple.
****************************************************************** */
template <class FunctionBase>
void PK_DomainFunctionSimple<FunctionBase>::Compute(double t0, double t1)
{
  if (unique_specs_.size() == 0) return;

  // create the input tuple (time + space)
  int dim = mesh_->space_dimension();
  std::vector<double> args(1 + dim);

  for (auto uspec = unique_specs_.at(kind_)->begin(); uspec != unique_specs_.at(kind_)->end(); ++uspec) {
    args[0] = t1;
    Teuchos::RCP<MeshIDs> ids = (*uspec)->second;
    // uspec->first is a RCP<Spec>, Spec's second is an RCP to the function.
    int nfun = (*uspec)->first->second->size();
    std::vector<double> val_vec(nfun);  

    for (auto c = ids->begin(); c != ids->end(); ++c) {
      auto xc = PKUtils_EntityCoordinates(*c, kind_, *mesh_);
      for (int i = 0; i != dim; ++i) args[i + 1] = xc[i];
      
      for (int i = 0; i < nfun; ++i) {
        val_vec[i] = (*(*uspec)->first->second)(args)[i];
      }
      value_[*c] = val_vec;
    }

    if (submodel_ == "integrated source") {
      double dt = t1 - t0;
 
      args[0] = t0;
      for (auto c = ids->begin(); c != ids->end(); ++c) {
        auto xc = PKUtils_EntityCoordinates(*c, kind_, *mesh_);

        for (int i = 0; i != dim; ++i) args[i + 1] = xc[i];
        for (int i = 0; i < nfun; ++i) {
          value_[*c][i] -= (*(*uspec)->first->second)(args)[i];
          if (dt > 0.0) value_[*c][i] /= dt;
        }
      }
    }
  }
}

}  // namespace Amanzi

#endif
