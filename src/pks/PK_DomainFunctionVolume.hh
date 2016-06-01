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
  PK_DomainFunctionVolume(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) :
      UniqueMeshFunction(mesh) {};
  ~PK_DomainFunctionVolume() {};

  // member functions
  void Init(const Teuchos::ParameterList& plist);

  // required member functions
  virtual void Compute(double t0, double t1);
  virtual std::string name() { return "volume"; }

 protected:
  using FunctionBase::value_;
  using FunctionBase::domain_volume_;

 private:
  std::string submodel_;
  std::map<AmanziMesh::Entity_kind, std::vector<double> > measuare_;
};


/* ******************************************************************
* Initialization adds a single function to the list of unique specs.
****************************************************************** */
template <class FunctionBase>
void PK_DomainFunctionVolume<FunctionBase>::Init(const Teuchos::ParameterList& plist)
{
  submodel_ = "rate";
  if (plist.isParameter("submodel"))
    submodel_ = plist.get<std::string>("submodel");

  std::vector<std::string> regions = plist.get<Teuchos::Array<std::string> >("regions").toVector();

  Teuchos::RCP<Amanzi::MultiFunction> f;
  try {
    Teuchos::ParameterList flist = plist.sublist("sink");
    f = Teuchos::rcp(new MultiFunction(flist));
  } catch (Errors::Message& msg) {
    Errors::Message m;
    m << "error in source sublist : " << msg.what();
    Exceptions::amanzi_throw(m);
  }

  // Add this source specification to the domain function.
  Teuchos::RCP<Domain> domain = Teuchos::rcp(new Domain(regions, AmanziMesh::CELL));
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

  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);

  for (UniqueSpecList::const_iterator uspec = unique_specs_[AmanziMesh::CELL]->begin();
       uspec != unique_specs_[AmanziMesh::CELL]->end(); ++uspec) {

    Teuchos::RCP<MeshIDs> ids = (*uspec)->second;

    // calculate physical volume of region defined by domain. 
    domain_volume_ = 0.0;
    for (MeshIDs::const_iterator c = ids->begin(); c != ids->end(); ++c) {
      if (*c < ncells_owned) domain_volume_ += mesh_->cell_volume(*c);
    }
    double tmp(domain_volume_);
    mesh_->get_comm()->SumAll(&tmp, &domain_volume_, 1);

    args[0] = t1;
    for (MeshIDs::const_iterator c = ids->begin(); c != ids->end(); ++c) {
      const AmanziGeometry::Point& xc = mesh_->cell_centroid(*c);
      for (int i = 0; i != dim; ++i) args[i + 1] = xc[i];

      // uspec->first is a RCP<Spec>, Spec's second is an RCP to the function.
      value_[*c] = (*(*uspec)->first->second)(args)[0] / domain_volume_;
    }

    if (submodel_ == "integrated source") {
      double dt = t1 - t0;
      if (dt > 0.0) dt = 1.0 / dt;

      args[0] = t0;
      for (MeshIDs::const_iterator c = ids->begin(); c != ids->end(); ++c) {
        const AmanziGeometry::Point& xc = mesh_->cell_centroid(*c);
        for (int i = 0; i != dim; ++i) args[i + 1] = xc[i];

        value_[*c] -= (*(*uspec)->first->second)(args)[0] / domain_volume_;
        value_[*c] *= dt;
      }
    }
  }
}

}  // namespace Amanzi

#endif

