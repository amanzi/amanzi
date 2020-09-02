/*
  Process Kernels 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Scenario 1 (data distribution model is "volume")
  The total source Q is given on input. The uniform source 
  distribution is employed. The source density in cell c
  calculated as (Q / V_D) * vol_fraction(c) 

  Scenario 2 (data distribution model is "none" or "simple")
  The local source function Q(x) is given on input. The source 
  source density in cell c is calculate as Q(x_c) * vol_fraction(c)
  where x_c is the centroid of volume_fraction.
*/

#ifndef AMANZI_PK_DOMAIN_FUNCTION_VOLUME_FRACTION_HH_
#define AMANZI_PK_DOMAIN_FUNCTION_VOLUME_FRACTION_HH_

#include <string>
#include <vector>

#include "Teuchos_RCP.hpp"

#include "CommonDefs.hh"
#include "MaterialMeshFunction.hh"
#include "Mesh.hh"
#include "PK_DomainFunction.hh"

namespace Amanzi {

template <class FunctionBase>
class PK_DomainFunctionVolumeFraction : public FunctionBase,
                                        public Functions::MaterialMeshFunction {
 public:
  PK_DomainFunctionVolumeFraction(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                                  AmanziMesh::Entity_kind kind) :
      MaterialMeshFunction(mesh),
      kind_(kind) {};
  ~PK_DomainFunctionVolumeFraction() {};

  // member functions
  void Init(const Teuchos::ParameterList& plist, const std::string& keyword);

  // required member functions
  virtual void Compute(double t0, double t1);
  virtual std::string name() const { return "volume fraction"; }

 protected:
  using FunctionBase::value_;
  using FunctionBase::domain_volume_;
  using FunctionBase::keyword_;

 private:
  std::string model_, submodel_;
  AmanziMesh::Entity_kind kind_;
};


/* ******************************************************************
* Initialization adds a single function to the list of specs.
****************************************************************** */
template <class FunctionBase>
void PK_DomainFunctionVolumeFraction<FunctionBase>::Init(
    const Teuchos::ParameterList& plist, const std::string& keyword)
{
  keyword_ = keyword;

  model_ = plist.get<std::string>("spatial distribution method");

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
void PK_DomainFunctionVolumeFraction<FunctionBase>::Compute(double t0, double t1)
{
  // create the input tuple (time + space)
  int dim = (*mesh_).space_dimension();
  std::vector<double> args(1 + dim);

  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);

  for (MaterialSpecList::const_iterator mspec = material_specs_.at(kind_)->begin();
       mspec != material_specs_.at(kind_)->end(); ++mspec) {

    Teuchos::RCP<MaterialMesh> ids = (*mspec)->second;

    // calculate physical volume of region.
    if (model_ == "volume") {
      domain_volume_ = 0.0;
      for (MaterialMesh::const_iterator it = ids->begin(); it != ids->end(); ++it) {
        int c = it->first;
        double vofs = it->second;
        if (c < ncells_owned) domain_volume_ += 
          ((kind_ == AmanziMesh::CELL) ? mesh_->cell_volume(c) : mesh_->face_area(c)) * vofs;
      }
      double tmp = domain_volume_;
      mesh_->get_comm()->SumAll(&tmp, &domain_volume_, 1);
    } else {
      domain_volume_ = 1.0;
    }

    int nfun = (*mspec)->first->second->size();
    std::vector<double> val_vec(nfun); 

    args[0] = t1;
    for (MaterialMesh::const_iterator it = ids->begin(); it != ids->end(); ++it) {
      int c = it->first;
      double vofs = it->second;

      auto xc = PKUtils_EntityCoordinates(c, kind_, *mesh_);
      for (int i = 0; i != dim; ++i) args[i + 1] = xc[i];

      // mspec->first is a RCP<Spec>, Spec's second is an RCP to the function.
      for (int i = 0; i < nfun; ++i) {
        val_vec[i] = (*(*mspec)->first->second)(args)[i] * vofs / domain_volume_;
      }
      value_[c] = val_vec;
    }

    if (submodel_ == "integrated source") {
      double dt = t1 - t0;

      args[0] = t0;
      for (MaterialMesh::const_iterator it = ids->begin(); it != ids->end(); ++it) {
        int c = it->first;
        double vofs = it->second;

        auto xc = PKUtils_EntityCoordinates(c, kind_, *mesh_);
        for (int i = 0; i != dim; ++i) args[i + 1] = xc[i];

        for (int i = 0; i < nfun; ++i) {
          value_[c][i] -= (*(*mspec)->first->second)(args)[i] * vofs / domain_volume_;
          if (dt > 0.0) value_[c][i] /= dt;
        }
      }
    }
  }
}

}  // namespace Amanzi

#endif

