/*
  Process Kernels 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
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

 private:
  std::string submodel_;
};


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

    double domain_volume = 0.0;
    Teuchos::RCP<MeshIDs> ids = (*uspec)->second;

    // calculate physical volume of region.
    for (MeshIDs::const_iterator c = ids->begin(); c != ids->end(); ++c) {
      if (*c < ncells_owned) domain_volume += mesh_->cell_volume(*c);
    }
    double volume_tmp = domain_volume;
    mesh_->get_comm()->SumAll(&volume_tmp, &domain_volume, 1);

    args[0] = t1;
    for (MeshIDs::const_iterator c = ids->begin(); c != ids->end(); ++c) {
      const AmanziGeometry::Point& xc = mesh_->cell_centroid(*c);
      for (int i = 0; i != dim; ++i) args[i + 1] = xc[i];

      // uspec->first is a RCP<Spec>, Spec's second is an RCP to the function.
      value_[*c] = (*(*uspec)->first->second)(args)[0] / domain_volume;
    }

    if (submodel_ == "integrated source") {
      double dt = t1 - t0;
      if (dt > 0.0) dt = 1.0 / dt;

      args[0] = t0;
      for (MeshIDs::const_iterator c = ids->begin(); c != ids->end(); ++c) {
        const AmanziGeometry::Point& xc = mesh_->cell_centroid(*c);
        for (int i = 0; i != dim; ++i) args[i + 1] = xc[i];

        value_[*c] -= (*(*uspec)->first->second)(args)[0] / domain_volume;
        value_[*c] *= dt;
      }
    }
  }
}

}  // namespace Amanzi

#endif

