/*
  Process Kernels 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  The total source Q is given for each domain. The uniform source 
  distribution model is employed. The cell-based source density is 
  calculated as (Q / V_D) * vol_fraction. 
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
  PK_DomainFunctionVolumeFraction(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) :
      MaterialMeshFunction(mesh) {};
  ~PK_DomainFunctionVolumeFraction() {};

  // member functions
  void Init(const Teuchos::ParameterList& plist);

  // required member functions
  virtual void Compute(double t0, double t1);
  virtual std::string name() { return "volume fraction"; }

 protected:
  using FunctionBase::value_;
  using FunctionBase::domain_volume_;

 private:
  std::string submodel_;
};


/* ******************************************************************
* Initialization adds a single function to the list of specs.
****************************************************************** */
template <class FunctionBase>
void PK_DomainFunctionVolumeFraction<FunctionBase>::Init(const Teuchos::ParameterList& plist)
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
void PK_DomainFunctionVolumeFraction<FunctionBase>::Compute(double t0, double t1)
{
   // create the input tuple (time + space)
  int dim = (*mesh_).space_dimension();
  std::vector<double> args(1 + dim);

  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::OWNED);

  for (MaterialSpecList::const_iterator mspec = material_specs_[AmanziMesh::CELL]->begin();
       mspec != material_specs_[AmanziMesh::CELL]->end(); ++mspec) {

    Teuchos::RCP<MaterialMesh> ids = (*mspec)->second;

    // calculate physical volume of region.
    domain_volume_ = 0.0;
    for (MaterialMesh::const_iterator it = ids->begin(); it != ids->end(); ++it) {
      int c = it->first;
      double vofs = it->second;
      if (c < ncells_owned) domain_volume_ += mesh_->cell_volume(c) * vofs;
    }
    double tmp = domain_volume_;
    mesh_->get_comm()->SumAll(&tmp, &domain_volume_, 1);

    args[0] = t1;
    for (MaterialMesh::const_iterator it = ids->begin(); it != ids->end(); ++it) {
      int c = it->first;
      double vofs = it->second;

      const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);
      for (int i = 0; i != dim; ++i) args[i + 1] = xc[i];

      // mspec->first is a RCP<Spec>, Spec's second is an RCP to the function.
      // std::cout << xc << " " << vofs << std::endl;
      value_[c] = (*(*mspec)->first->second)(args)[0] * vofs / domain_volume_;
    }

    if (submodel_ == "integrated source") {
      double dt = t1 - t0;
      if (dt > 0.0) dt = 1.0 / dt;

      args[0] = t0;
      for (MaterialMesh::const_iterator it = ids->begin(); it != ids->end(); ++it) {
        int c = it->first;
        double vofs = it->second;

        const AmanziGeometry::Point& xc = mesh_->cell_centroid(c);
        for (int i = 0; i != dim; ++i) args[i + 1] = xc[i];

        value_[c] -= (*(*mspec)->first->second)(args)[0] * vofs / domain_volume_;
        value_[c] *= dt;
      }
    }
  }
}

}  // namespace Amanzi

#endif

