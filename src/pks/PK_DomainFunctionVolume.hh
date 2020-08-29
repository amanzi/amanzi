/*
  Process Kernels 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  The total source Q is given for each domain. The uniform source 
  distribution model is employed. The cell-based source density is 
  calculated as (Q / V_D). 
*/

#ifndef AMANZI_PK_DOMAIN_FUNCTION_VOLUME_HH_
#define AMANZI_PK_DOMAIN_FUNCTION_VOLUME_HH_

#include <string>
#include <vector>

#include "Teuchos_RCP.hpp"

#include "CommonDefs.hh"
#include "Mesh.hh"
#include "PK_DomainFunction.hh"

namespace Amanzi {

template <class FunctionBase>
class PK_DomainFunctionVolume : public FunctionBase,
                                public Functions::UniqueMeshFunction {
 public:
  PK_DomainFunctionVolume(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                          AmanziMesh::Entity_kind kind) :
      UniqueMeshFunction(mesh),
      kind_(kind) {};

  PK_DomainFunctionVolume(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                          const Teuchos::ParameterList& plist,
                          AmanziMesh::Entity_kind kind) :
      UniqueMeshFunction(mesh),
      kind_(kind) {};

  ~PK_DomainFunctionVolume() {};

  // member functions
  void Init(const Teuchos::ParameterList& plist, const std::string& keyword);

  // required member functions
  virtual void Compute(double t0, double t1);
  virtual std::string name() const { return "volume"; }

 protected:
  using FunctionBase::value_;
  using FunctionBase::domain_volume_;
  using FunctionBase::keyword_;

 private:
  std::string submodel_;
  AmanziMesh::Entity_kind kind_;
  std::map<AmanziMesh::Entity_kind, std::vector<double> > measure_;
};


/* ******************************************************************
* Initialization adds a single function to the list of unique specs.
****************************************************************** */
template <class FunctionBase>
void PK_DomainFunctionVolume<FunctionBase>::Init(
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
* Compute and distribute the result by volume.
****************************************************************** */
template <class FunctionBase>
void PK_DomainFunctionVolume<FunctionBase>::Compute(double t0, double t1)
{
   // create the input tuple (time + space)
  int dim = (*mesh_).space_dimension();
  std::vector<double> args(1 + dim);

  int nowned = mesh_->num_entities(kind_, AmanziMesh::Parallel_type::OWNED);

  for (auto uspec = unique_specs_.at(kind_)->begin(); uspec != unique_specs_.at(kind_)->end(); ++uspec) {
    Teuchos::RCP<MeshIDs> ids = (*uspec)->second;

    // calculate physical volume of region defined by domain. 
    domain_volume_ = 0.0;
    for (MeshIDs::const_iterator c = ids->begin(); c != ids->end(); ++c) {
      if (*c < nowned) domain_volume_ += 
          (kind_ == AmanziMesh::CELL) ? mesh_->cell_volume(*c) : mesh_->face_area(*c);
    }
    double tmp(domain_volume_);
    mesh_->get_comm()->SumAll(&tmp, &domain_volume_, 1);
    int nfun = (*uspec)->first->second->size();
    std::vector<double> val_vec(nfun);  

    args[0] = t1;
    for (auto c = ids->begin(); c != ids->end(); ++c) {
      auto xc = PKUtils_EntityCoordinates(*c, kind_, *mesh_);

      for (int i = 0; i != dim; ++i) args[i + 1] = xc[i];

      // uspec->first is a RCP<Spec>, Spec's second is an RCP to the function.
      for (int i = 0; i < nfun; ++i) {
        val_vec[i] = (*(*uspec)->first->second)(args)[i] / domain_volume_;
      }
      value_[*c] = val_vec;
    }

    if (submodel_ == "integrated source") {
      double dt = t1 - t0;
      if (dt > 0.0) dt = 1.0 / dt;

      args[0] = t0;
      for (auto c = ids->begin(); c != ids->end(); ++c) {
        auto xc = PKUtils_EntityCoordinates(*c, kind_, *mesh_);
        for (int i = 0; i != dim; ++i) args[i + 1] = xc[i];

        for (int i = 0; i < nfun; ++i) {
          value_[*c][i] -= (*(*uspec)->first->second)(args)[i] / domain_volume_;
          value_[*c][i] *= dt;
        }
      }
    }
  }
}

}  // namespace Amanzi

#endif

