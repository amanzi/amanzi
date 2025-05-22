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

#ifndef AMANZI_PK_DOMAIN_FUNCTION_SUBGRID_RETURN_HH_
#define AMANZI_PK_DOMAIN_FUNCTION_SUBGRID_RETURN_HH_

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

template <class FunctionBase>
class PK_DomainFunctionSubgridReturn : public FunctionBase, public Functions::UniqueMeshFunction {
 public:
  PK_DomainFunctionSubgridReturn(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                                 const Teuchos::ParameterList& plist)
    : FunctionBase(plist), UniqueMeshFunction(mesh, AmanziMesh::Parallel_kind::OWNED){};
  virtual ~PK_DomainFunctionSubgridReturn() = default;

  // member functions
  void Init(const Teuchos::ParameterList& plist, const std::string& keyword);

  // required member functions
  virtual void Compute(double t0, double t1) override;
  virtual DomainFunction_kind getType() const override
  {
    return DomainFunction_kind::SUBGRID_RETURN;
  }
  virtual void set_state(const Teuchos::RCP<State>& S) final { S_ = S; }

 protected:
  using FunctionBase::value_;
  using FunctionBase::keyword_;
  using Functions::UniqueMeshFunction::mesh_;
  Teuchos::RCP<State> S_;

 private:
  Key exchanged_suffix_, dset_;
  Tag exchanged_tag_;

  Key lwc_key_;
  Tag lwc_tag_;
};


/* ******************************************************************
* Initialization
****************************************************************** */
template <class FunctionBase>
void
PK_DomainFunctionSubgridReturn<FunctionBase>::Init(const Teuchos::ParameterList& plist,
                                                   const std::string& keyword)
{
  keyword_ = keyword;

  // get and check the model parameters
  Teuchos::ParameterList blist = plist.sublist("source function");
  Key domain_name = Keys::readDomain(blist, "surface", "surface");

  exchanged_suffix_ = blist.get<std::string>("subgrid field suffix");
  exchanged_tag_ = Keys::readTag(blist, "subgrid field");

  // domain set of the collection of subgrids
  dset_ = blist.get<std::string>("subgrid domain set", "subgrid");

  lwc_key_ = Keys::readKey(blist, domain_name, "liquid water content", "water_content");
  lwc_tag_ = Keys::readTag(blist, "liquid water content");

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
  auto domain = Teuchos::rcp(new Domain(regions, AmanziMesh::Entity_kind::CELL));
  AddSpec(Teuchos::rcp(new Spec(domain, f)));
}


/* ******************************************************************
* Compute and distribute the result by SubgridReturn.
****************************************************************** */
template <class FunctionBase>
void
PK_DomainFunctionSubgridReturn<FunctionBase>::Compute(double t0, double t1)
{
  if (unique_specs_.size() == 0) return;

  // create the input tuple (time + space)
  int dim = mesh_->getSpaceDimension();
  std::vector<double> args(1 + dim);
  args[0] = t1;

  // get the map to convert to subgrid GID
  auto& map = mesh_->getMap(AmanziMesh::Entity_kind::CELL, false);

  S_->GetEvaluator(lwc_key_, lwc_tag_).Update(*S_, lwc_key_+"_subgrid_return");
  const auto& lwc = *S_->Get<CompositeVector>(lwc_key_, lwc_tag_).ViewComponent("cell");

  for (const auto& uspec : *unique_specs_.at(AmanziMesh::Entity_kind::CELL)) {
    Teuchos::RCP<MeshIDs> ids = uspec->second;

    // uspec->first is a RCP<Spec>, Spec's second is an RCP to the function.
    int nfun = uspec->first->second->size();

    // loop over all entities in the spec
    for (const AmanziMesh::Entity_ID& c : *ids) {
      const AmanziGeometry::Point& xc = mesh_->getCellCentroid(c);

      // set the xyz coords
      for (int i = 0; i != dim; ++i) args[i + 1] = xc[i];

      // uspec->first is a RCP<Spec>, Spec's second is an RCP to the function.
      std::vector<double> alpha(nfun);
      for (int i = 0; i < nfun; ++i) { alpha[i] = (*uspec->first->second)(args)[i]; }

      // find the subgrid gid to be integrated
      auto gid = map.GID(c);
      Key domain = Keys::getDomainInSet(dset_, gid);
      Key var = Keys::getKey(domain, exchanged_suffix_);

      // get the vector to be integrated
      const auto& vec_out = S_->Get<CompositeVector>(var, exchanged_tag_);
      const auto& vec_c = *vec_out.ViewComponent("cell", true);

      std::vector<double> val(nfun, 0.);

      // DO THE INTEGRAL: currently omega_i = 1/cv_sg?
      int ncells_sg = vec_out.Mesh()->getNumEntities(AmanziMesh::Entity_kind::CELL,
                                                     AmanziMesh::Parallel_kind::ALL);
      for (int c_sg = 0; c_sg != ncells_sg; ++c_sg) {
        for (int k = 0; k != nfun; ++k) { val[k] += vec_c[k][c_sg] * alpha[k]; }
      }
      for (int k = 0; k != nfun; ++k) {
        val[k] *= lwc[0][c] / mesh_->getCellVolume(c) / ncells_sg;
      }
      value_[c] = val;
    }
  }
}

} // namespace Amanzi

#endif
