/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*!

To be written.

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

template<class FunctionBase>
class PK_DomainFunctionSimple
  : public FunctionBase
  , public Functions::UniqueMeshFunction {
 public:
  PK_DomainFunctionSimple(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                          AmanziMesh::Entity_kind entity_kind,
                          bool ghosted)
    : UniqueMeshFunction(mesh,
                         ghosted ? AmanziMesh::Parallel_kind::ALL :
                                   AmanziMesh::Parallel_kind::OWNED),
      kind_(entity_kind) {};

  PK_DomainFunctionSimple(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                          const Teuchos::ParameterList& plist,
                          AmanziMesh::Entity_kind entity_kind,
                          bool ghosted)
    : FunctionBase(plist),
      UniqueMeshFunction(mesh,
                         ghosted ? AmanziMesh::Parallel_kind::ALL :
                                   AmanziMesh::Parallel_kind::OWNED),
      kind_(entity_kind) {};

  ~PK_DomainFunctionSimple() {};

  // member functions
  void Init(const Teuchos::ParameterList& plist, const std::string& keyword);

  // required member functions
  virtual void Compute(double t0, double t1) override;
  virtual DomainFunction_kind getType() const override { return DomainFunction_kind::SIMPLE; }

 protected:
  using FunctionBase::keyword_;
  using FunctionBase::value_;

 private:
  std::string submodel_;
  AmanziMesh::Entity_kind kind_;
};


/* ******************************************************************
* Initialization adds a single function to the list of unique specs.
****************************************************************** */
template<class FunctionBase>
void
PK_DomainFunctionSimple<FunctionBase>::Init(const Teuchos::ParameterList& plist,
                                            const std::string& keyword)
{
  keyword_ = keyword;

  submodel_ = "rate";
  if (plist.isParameter("submodel") ) submodel_ = plist.get<std::string>("submodel");
  std::vector<std::string> regions = plist.get<Teuchos::Array<std::string>>("regions").toVector();

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
template<class FunctionBase>
void
PK_DomainFunctionSimple<FunctionBase>::Compute(double t0, double t1)
{
  if (unique_specs_.size() == 0) return;

  // create the input tuple (time + space)
  double dt = t1 - t0;
  int dim = mesh_->getSpaceDimension();
  std::vector<double> args(1 + dim);

  for (auto uspec = unique_specs_.at(kind_) ->begin(); uspec != unique_specs_.at(kind_)->end();
       ++uspec) {
    Teuchos::RCP<MeshIDs> ids = (*uspec)->second;
    // uspec->first is a RCP<Spec>, Spec's second is an RCP to the function.
    int nfun = (*uspec)->first->second->size();

    for (auto c = ids->begin(); c != ids->end(); ++c) {
      value_[*c].resize(nfun);

      args[0] = t1;
      auto xc = PKUtils_EntityCoordinates(*c, kind_, *mesh_);
      for (int i = 0; i != dim; ++i) args[i + 1] = xc[i];

      const double* val1 = (*(*uspec)->first->second)(args);
      for (int i = 0; i < nfun; ++i) value_[*c][i] = val1[i];

      if (submodel_ == "integrated source") {
        args[0] = t0;
        const double* val0 = (*(*uspec)->first->second)(args);
        for (int i = 0; i < nfun; ++i) {
          value_[*c][i] -= val0[i];
          if (dt > 0.0) value_[*c][i] /= dt;
        }
      }
    }
  }
}

} // namespace Amanzi

#endif
