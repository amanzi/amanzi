/*
  Process Kernels 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (coonet@ornl.gov)
*/

#ifndef AMANZI_PK_DOMAIN_FUNCTION_FIRSTORDER_EXCHANGE_HH_
#define AMANZI_PK_DOMAIN_FUNCTION_FIRSTORDER_EXCHANGE_HH_

#include <string>
#include <vector>

#include "Epetra_Vector.h"
#include "Teuchos_RCP.hpp"

#include "CommonDefs.hh"
#include "DenseVector.hh"
#include "Key.hh"
#include "Mesh.hh"
#include "Tag.hh"
#include "UniqueMeshFunction.hh"


namespace Amanzi {

template <class FunctionBase>
class PK_DomainFunctionFirstOrderExchange : public FunctionBase,
                                            public Functions::UniqueMeshFunction {
 public:
  PK_DomainFunctionFirstOrderExchange(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                                      const Teuchos::ParameterList& plist,
                                      AmanziMesh::Entity_kind kind)
    : FunctionBase(plist), UniqueMeshFunction(mesh), kind_(kind){};
  virtual ~PK_DomainFunctionFirstOrderExchange() = default;

  // member functions
  void Init(const Teuchos::ParameterList& plist, const std::string& keyword);

  virtual void set_state(const Teuchos::RCP<State>& S) { S_ = S; }

  // required member functions
  virtual void Compute(double t0, double t1);
  virtual std::string name() const { return "first order exchange"; }

 protected:
  using FunctionBase::value_;
  using FunctionBase::keyword_;

  Teuchos::RCP<const State> S_;

 private:
  AmanziMesh::Entity_kind kind_;
  Key tcc_key_;
  Tag tcc_copy_;

  Key saturation_key_, porosity_key_, molar_density_key_;
  Tag saturation_copy_, porosity_copy_, molar_density_copy_;
};


/* ******************************************************************
* Initialization adds a single function to the list of unique specs.
****************************************************************** */
template <class FunctionBase>
void
PK_DomainFunctionFirstOrderExchange<FunctionBase>::Init(const Teuchos::ParameterList& plist,
                                                        const std::string& keyword)
{
  keyword_ = keyword;

  // get the model parameters
  Teuchos::ParameterList blist = plist.sublist("source function");
  std::string domain_name = blist.get<std::string>("domain name", "domain");

  tcc_key_ = Keys::readKey(
    blist, domain_name, "total component concentration", "total_component_concentration");
  tcc_copy_ = Keys::readTag(blist, "total component concentration copy");

  saturation_key_ = Keys::readKey(blist, domain_name, "saturation liquid", "ponded_depth");
  saturation_copy_ = Keys::readTag(blist, "saturation liquid copy");

  porosity_key_ = Keys::readKey(blist, domain_name, "porosity", "porosity");
  porosity_copy_ = Keys::readTag(blist, "porosity copy");

  molar_density_key_ = Keys::readKey(blist, domain_name, "molar density", "molar_density_liquid");
  molar_density_copy_ = Keys::readTag(blist, "molar density copy");

  // get and check the regions
  auto regions = plist.get<Teuchos::Array<std::string>>("regions").toVector();

  // get the function for alpha
  Teuchos::RCP<Amanzi::MultiFunction> f;
  try {
    Teuchos::ParameterList flist = blist.sublist("function");
    f = Teuchos::rcp(new MultiFunction(flist));
  } catch (Errors::Message& msg) {
    Errors::Message m;
    m << "error in source sublist : " << msg.what();
    Exceptions::amanzi_throw(m);
  }

  // Add this source specification to the domain function.
  auto domain = Teuchos::rcp(new Domain(regions, kind_));
  AddSpec(Teuchos::rcp(new Spec(domain, f)));
}


/* ******************************************************************
* Compute and distribute the result by FirstOrderExchange.
****************************************************************** */
template <class FunctionBase>
void
PK_DomainFunctionFirstOrderExchange<FunctionBase>::Compute(double t0, double t1)
{
  if (unique_specs_.size() == 0) return;

  // create the input tuple (time + space)
  int dim = mesh_->space_dimension();
  std::vector<double> args(1 + dim);

  // get the tcc vector
  const auto& tcc = *S_->Get<CompositeVector>(tcc_key_, tcc_copy_).ViewComponent("cell");
  const auto& ws_ =
    *S_->Get<CompositeVector>(saturation_key_, saturation_copy_).ViewComponent("cell");
  const auto& phi_ = *S_->Get<CompositeVector>(porosity_key_, porosity_copy_).ViewComponent("cell");
  const auto& mol_dens_ =
    *S_->Get<CompositeVector>(molar_density_key_, molar_density_copy_).ViewComponent("cell");

  for (UniqueSpecList::const_iterator uspec = unique_specs_.at(kind_)->begin();
       uspec != unique_specs_.at(kind_)->end();
       ++uspec) {
    args[0] = t1;
    Teuchos::RCP<MeshIDs> ids = (*uspec)->second;
    // uspec->first is a RCP<Spec>, Spec's second is an RCP to the function.
    int nfun = (*uspec)->first->second->size();
    std::vector<double> val_vec(nfun, 0.);

    for (auto c = ids->begin(); c != ids->end(); ++c) {
      auto xc = PKUtils_EntityCoordinates(*c, kind_, *mesh_);

      for (int i = 0; i != dim; ++i) args[i + 1] = xc[i];

      // uspec->first is a RCP<Spec>, Spec's second is an RCP to the function.
      for (int i = 0; i < nfun; ++i) {
        val_vec[i] = -(*(*uspec)->first->second)(args)[i] * tcc[i][*c] * ws_[0][*c] * phi_[0][*c] *
                     mol_dens_[0][*c];
      }
      value_[*c] = val_vec;
    }
  }
}

} // namespace Amanzi

#endif
