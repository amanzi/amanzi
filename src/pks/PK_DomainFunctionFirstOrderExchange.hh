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

template<class FunctionBase>
class PK_DomainFunctionFirstOrderExchange
  : public FunctionBase
  , public Functions::UniqueMeshFunction {
 public:
  PK_DomainFunctionFirstOrderExchange(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                                      const Teuchos::ParameterList& plist,
                                      AmanziMesh::Entity_kind kind)
    : FunctionBase(plist),
      UniqueMeshFunction(mesh, AmanziMesh::Parallel_kind::OWNED),
      kind_(kind) {};
  virtual ~PK_DomainFunctionFirstOrderExchange() = default;

  // member functions
  void Init(const Teuchos::ParameterList& plist, const std::string& keyword);

  virtual void set_state(const Teuchos::RCP<State>& S) final { S_ = S; }

  // required member functions
  virtual void Compute(double t0, double t1) override;
  virtual DomainFunction_kind getType() const override
  {
    return DomainFunction_kind::FIRST_ORDER_EXCHANGE;
  }

 protected:
  using FunctionBase::keyword_;
  using FunctionBase::value_;

  Teuchos::RCP<State> S_;

 private:
  AmanziMesh::Entity_kind kind_;
  Key exchanged_key_;
  Tag exchanged_tag_;

  Key lwc_key_;
  Tag lwc_tag_;
};


/* ******************************************************************
* Initialization adds a single function to the list of unique specs.
****************************************************************** */
template<class FunctionBase>
void
PK_DomainFunctionFirstOrderExchange<FunctionBase>::Init(const Teuchos::ParameterList& plist,
                                                        const std::string& keyword)
{
  Teuchos::ParameterList blist = plist.sublist("source function");

  Key domain_name = Keys::readDomain(blist, "domain", "domain");
  exchanged_key_ = Keys::readKey(blist, domain_name, "exchanged quantity", "mole_fraction");
  exchanged_tag_ = Keys::readTag(blist, "exchanged quantity");

  lwc_key_ = Keys::readKey(blist, domain_name, "liquid water content", "water_content");
  lwc_tag_ = Keys::readTag(blist, "liquid water content", exchanged_tag_);

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
template<class FunctionBase>
void
PK_DomainFunctionFirstOrderExchange<FunctionBase>::Compute(double t0, double t1)
{
  std::cout << "FirstOrderExchange" << std::endl;
  if (unique_specs_.size() == 0) return;

  // create the input tuple (time + space)
  int dim = mesh_->getSpaceDimension();
  std::vector<double> args(1 + dim);
  args[0] = t1;

  // get the mole fraction, water contents
  S_->GetEvaluator(exchanged_key_, exchanged_tag_)
    .Update(*S_, exchanged_key_ + "_first_order_exchange");
  S_->GetEvaluator(lwc_key_, lwc_tag_).Update(*S_, lwc_key_ + "_first_order_exchange");
  const auto& ex = *S_->Get<CompositeVector>(exchanged_key_, exchanged_tag_).ViewComponent("cell");
  const auto& lwc = *S_->Get<CompositeVector>(lwc_key_, lwc_tag_).ViewComponent("cell");

  for (const auto& uspec : *unique_specs_.at(kind_)) {
    // uspec is a RCP<pair<RCP<Spec>>, RCP<MeshIDs>>
    Teuchos::RCP<MeshIDs> ids = uspec->second;
    int nfun = uspec->first->second->size();
    std::vector<double> val_vec(nfun, 0.);

    for (auto c : *ids) {
      const auto& xc = mesh_->getCentroid(kind_, c);
      for (int i = 0; i != dim; ++i) args[i + 1] = xc[i];

      // uspec->first is a RCP<Spec>, Spec's second is an RCP to the function.
      for (int i = 0; i < nfun; ++i) {
        val_vec[i] =
          -(*uspec->first->second)(args)[i] * ex[i][c] * lwc[0][c] / mesh_->getCellVolume(c);
      }
      value_[c] = val_vec;
    }
  }
}

} // namespace Amanzi

#endif
