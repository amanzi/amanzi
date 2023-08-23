/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatsky (dasvyat@lanl.gov)
*/

/*
  Process Kernels

*/

#ifndef AMANZI_PK_DOMAIN_FUNCTION_SIMPLEWELL_HH_
#define AMANZI_PK_DOMAIN_FUNCTION_SIMPLEWELL_HH_

#include <string>
#include <vector>

#include "Teuchos_RCP.hpp"

#include "CommonDefs.hh"
#include "Mesh.hh"
#include "State.hh"
#include "PK_DomainFunction.hh"

namespace Amanzi {

template <class FunctionBase>
class PK_DomainFunctionSimpleWell : public FunctionBase, public Functions::UniqueMeshFunction {
 public:
  PK_DomainFunctionSimpleWell(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh)
    : UniqueMeshFunction(mesh){};

  PK_DomainFunctionSimpleWell(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                              const Teuchos::ParameterList& plist)
    : UniqueMeshFunction(mesh){};

  ~PK_DomainFunctionSimpleWell(){};

  // member functions
  void Init(const Teuchos::ParameterList& plist,
            const std::string& keyword,
            const Teuchos::RCP<const State>& S);

  // required member functions
  virtual void Compute(double t0, double t1);
  virtual std::string name() const { return "simple well"; }

 protected:
  using FunctionBase::value_;
  using FunctionBase::domain_volume_;
  using FunctionBase::keyword_;

  Teuchos::RCP<const State> S_;

 private:
  std::string submodel_;
  AmanziMesh::Entity_kind kind_;
  std::map<AmanziMesh::Entity_kind, std::vector<double>> measure_;
  Teuchos::RCP<Domain> domain_;
  AmanziGeometry::Point gravity_;
  double depth_, rho_;
  Key well_index_key_;
};


/* ******************************************************************
* Initialization adds a single function to the list of unique specs.
****************************************************************** */
template <class FunctionBase>
void
PK_DomainFunctionSimpleWell<FunctionBase>::Init(const Teuchos::ParameterList& plist,
                                                const std::string& keyword,
                                                const Teuchos::RCP<const State>& S)
{
  keyword_ = keyword;
  S_ = S;
  kind_ = AmanziMesh::Entity_kind::CELL;
  Teuchos::ParameterList well_list = plist.sublist(keyword);
  submodel_ = "rate";
  if (well_list.isParameter("submodel")) submodel_ = well_list.get<std::string>("submodel");

  // well_index_key_ = plist.get<std::string>("well_index_key", "well_index");
  well_index_key_ = "well_index";

  std::vector<std::string> regions = plist.get<Teuchos::Array<std::string>>("regions").toVector();

  Teuchos::RCP<Amanzi::MultiFunction> f;
  Teuchos::ParameterList flist = well_list.sublist(submodel_);

  try {
    f = Teuchos::rcp(new MultiFunction(flist));
  } catch (Errors::Message& msg) {
    Errors::Message m;
    m << "error in source sublist : " << msg.what();
    Exceptions::amanzi_throw(m);
  }

  if (submodel_ == "bhp") {
    depth_ = well_list.get<double>("depth");
    gravity_ = S_->Get<AmanziGeometry::Point>("gravity", Tags::DEFAULT);
    rho_ = S->Get<double>("const_fluid_density", Tags::DEFAULT);
  }

  // Add this source specification to the domain function.
  domain_ = Teuchos::rcp(new Domain(regions, kind_));
  AddSpec(Teuchos::rcp(new Spec(domain_, f)));
}


/* ******************************************************************
* Compute and distribute the result by volume.
****************************************************************** */
template <class FunctionBase>
void
PK_DomainFunctionSimpleWell<FunctionBase>::Compute(double t0, double t1)
{
  // create the input tuple (time + space)
  int dim = mesh_->getSpaceDimension();
  std::vector<double> args(1 + dim);

  int nowned = mesh_->getNumEntities(kind_, AmanziMesh::Parallel_type::OWNED);

  if (submodel_ == "rate") {
    for (auto uspec = unique_specs_.at(kind_)->begin(); uspec != unique_specs_.at(kind_)->end();
         ++uspec) {
      Teuchos::RCP<MeshIDs> ids = (*uspec)->second;

      // calculate physical volume of region defined by domain.
      domain_volume_ = 0.0;
      for (MeshIDs::const_iterator c = ids->begin(); c != ids->end(); ++c) {
        if (*c < nowned)
          domain_volume_ += (kind_ == AmanziMesh::Entity_kind::CELL) ? mesh_->getCellVolume(*c) :
                                                                       mesh_->getFaceArea(*c);
      }
      double tmp(domain_volume_);
      mesh_->getComm()->SumAll(&tmp, &domain_volume_, 1);
      int nfun = (*uspec)->first->second->size();
      std::vector<double> val_vec(nfun);

      args[0] = t1;
      for (auto c = ids->begin(); c != ids->end(); ++c) {
        auto xc = PKUtils_EntityCoordinates(*c, kind_, *mesh_);
        for (int i = 0; i != dim; ++i) args[i + 1] = xc[i];

        // uspec->first is a RCP<Spec>, Spec's second is an RCP to the function.
        for (int i = 0; i < nfun; ++i)
          val_vec[i] = (*(*uspec)->first->second)(args)[i] / domain_volume_;
        value_[*c] = val_vec;
      }
    }

  } else if (submodel_ == "bhp") {
    double g = fabs(gravity_[dim - 1]);
    const Epetra_MultiVector& wi =
      *S_->Get<CompositeVector>("well_index", Tags::DEFAULT).ViewComponent("cell");

    for (auto uspec = unique_specs_.at(kind_)->begin(); uspec != unique_specs_.at(kind_)->end();
         ++uspec) {
      Teuchos::RCP<MeshIDs> ids = (*uspec)->second;

      int nfun = (*uspec)->first->second->size();
      std::vector<double> val_vec(nfun);
      args[0] = t1;
      for (auto c = ids->begin(); c != ids->end(); ++c) {
        auto xc = PKUtils_EntityCoordinates(*c, kind_, *mesh_);
        for (int i = 0; i != dim; ++i) args[i + 1] = xc[i];

        double bhp;
        for (int i = 0; i < nfun; ++i) {
          bhp = (*(*uspec)->first->second)(args)[i] + rho_ * g * (depth_ - xc[dim - 1]);
          val_vec[i] = bhp * wi[0][*c] / mesh_->getCellVolume(*c);
        }
        value_[*c] = val_vec;
      }
    }
  }
}

} // namespace Amanzi

#endif
