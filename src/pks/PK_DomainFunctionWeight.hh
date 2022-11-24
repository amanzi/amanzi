/*
  Process Kernels 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  The total source Q is given for each domain. A weighted source 
  distribution model is employed. The cell-based source density is 
  calculated as (Q / W_D) * weight, where W_D is the weighted 
  domain volume. The weight is defined globally, for the whole 
  computational domain.
*/

#ifndef AMANZI_PK_DOMAIN_FUNCTION_WEIGHT_HH_
#define AMANZI_PK_DOMAIN_FUNCTION_WEIGHT_HH_

#include <string>
#include <vector>

#include "Epetra_Vector.h"
#include "Teuchos_RCP.hpp"

#include "CommonDefs.hh"
#include "Mesh.hh"
#include "PK_DomainFunction.hh"

namespace Amanzi {

template <class FunctionBase>
class PK_DomainFunctionWeight : public FunctionBase, public Functions::UniqueMeshFunction {
 public:
  PK_DomainFunctionWeight(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh)
    : UniqueMeshFunction(mesh){};
  ~PK_DomainFunctionWeight(){};

  // member functions
  void Init(const Teuchos::ParameterList& plist,
            const std::string& keyword,
            Teuchos::RCP<const Epetra_Vector> weight);

  // required member functions
  virtual void Compute(double t0, double t1);
  virtual std::string name() const { return "weight"; }

 protected:
  using FunctionBase::value_;
  using FunctionBase::domain_volume_;
  using FunctionBase::keyword_;

 private:
  std::string submodel_;
  Teuchos::RCP<const Epetra_Vector> weight_;
};


/* ******************************************************************
* Initialization adds a single function to the list of unique specs.
****************************************************************** */
template <class FunctionBase>
void
PK_DomainFunctionWeight<FunctionBase>::Init(const Teuchos::ParameterList& plist,
                                            const std::string& keyword,
                                            Teuchos::RCP<const Epetra_Vector> weight)
{
  AMANZI_ASSERT(weight != Teuchos::null);

  keyword_ = keyword;

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
  Teuchos::RCP<Domain> domain = Teuchos::rcp(new Domain(regions, AmanziMesh::CELL));
  AddSpec(Teuchos::rcp(new Spec(domain, f)));

  weight_ = weight;
}


/* ******************************************************************
* Compute and distribute the result by volume.
****************************************************************** */
template <class FunctionBase>
void
PK_DomainFunctionWeight<FunctionBase>::Compute(double t0, double t1)
{
  double dt = t1 - t0;
  if (dt > 0.0) dt = 1.0 / dt;

  // create the input tuple (time + space)
  int dim = mesh_->space_dimension();
  std::vector<double> args(1 + dim);

  int ncells_owned = mesh_->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::OWNED);

  for (UniqueSpecList::const_iterator uspec = unique_specs_[AmanziMesh::CELL]->begin();
       uspec != unique_specs_[AmanziMesh::CELL]->end();
       ++uspec) {
    domain_volume_ = 0.0;
    double vol, weight_volume = 0.0;
    Teuchos::RCP<MeshIDs> ids = (*uspec)->second;

    for (MeshIDs::const_iterator c = ids->begin(); c != ids->end(); ++c) {
      if (*c < ncells_owned) {
        vol = mesh_->cell_volume(*c);
        domain_volume_ += vol;
        weight_volume += vol * (*weight_)[*c];
      }
    }
    double result[2], tmp[2] = { domain_volume_, weight_volume };
    mesh_->get_comm()->SumAll(tmp, result, 2);
    domain_volume_ = result[0];
    weight_volume = result[1];
    int nfun = (*uspec)->first->second->size();
    std::vector<double> val_vec(nfun);

    args[0] = t1;
    for (MeshIDs::const_iterator c = ids->begin(); c != ids->end(); ++c) {
      const AmanziGeometry::Point& xc = mesh_->cell_centroid(*c);
      for (int i = 0; i != dim; ++i) args[i + 1] = xc[i];
      for (int i = 0; i < nfun; ++i) {
        val_vec[i] = (*(*uspec)->first->second)(args)[i] * (*weight_)[*c] / weight_volume;
      }
      value_[*c] = val_vec;
    }

    if (submodel_ == "integrated source") {
      args[0] = t0;
      for (MeshIDs::const_iterator c = ids->begin(); c != ids->end(); ++c) {
        const AmanziGeometry::Point& xc = mesh_->cell_centroid(*c);
        for (int i = 0; i != dim; ++i) args[i + 1] = xc[i];
        for (int i = 0; i < nfun; ++i) {
          value_[*c][i] -= (*(*uspec)->first->second)(args)[i] * (*weight_)[*c] / weight_volume;
          value_[*c][i] *= dt;
        }
      }
    }
  }
}

} // namespace Amanzi

#endif
