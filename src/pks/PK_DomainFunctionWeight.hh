/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*!

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

#include "Epetra_MultiVector.h"
#include "Teuchos_RCP.hpp"

#include "CommonDefs.hh"
#include "Mesh.hh"
#include "PK_DomainFunction.hh"

namespace Amanzi {

template <class FunctionBase>
class PK_DomainFunctionWeight : public FunctionBase, public Functions::UniqueMeshFunction {
 public:
  PK_DomainFunctionWeight(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                          AmanziMesh::Entity_kind kind)
    : UniqueMeshFunction(mesh), kind_(kind){};
  ~PK_DomainFunctionWeight(){};

  // member functions
  void Init(const Teuchos::ParameterList& plist,
            const std::string& keyword,
            Teuchos::RCP<const Epetra_MultiVector> weight);

  // required member functions
  virtual void Compute(double t0, double t1);
  virtual std::string name() const { return "weight"; }

 protected:
  using FunctionBase::value_;
  using FunctionBase::domain_volume_;
  using FunctionBase::keyword_;

  Teuchos::RCP<const Epetra_MultiVector> weight_;

 private:
  std::string submodel_;
  AmanziMesh::Entity_kind kind_;
};


/* ******************************************************************
* Initialization adds a single function to the list of unique specs.
****************************************************************** */
template <class FunctionBase>
void
PK_DomainFunctionWeight<FunctionBase>::Init(const Teuchos::ParameterList& plist,
                                            const std::string& keyword,
                                            Teuchos::RCP<const Epetra_MultiVector> weight)
{
  AMANZI_ASSERT(weight != Teuchos::null);

  keyword_ = keyword;
  weight_ = weight;

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
  Teuchos::RCP<Domain> domain = Teuchos::rcp(new Domain(regions, AmanziMesh::Entity_kind::CELL));
  AddSpec(Teuchos::rcp(new Spec(domain, f)));
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
  int dim = mesh_->getSpaceDimension();
  std::vector<double> args(1 + dim);

  int ncells_owned =
    mesh_->getNumEntities(AmanziMesh::Entity_kind::CELL, AmanziMesh::Parallel_kind::OWNED);

  for (UniqueSpecList::const_iterator uspec = unique_specs_[AmanziMesh::Entity_kind::CELL]->begin();
       uspec != unique_specs_[AmanziMesh::Entity_kind::CELL]->end();
       ++uspec) {
    domain_volume_ = 0.0;
    double vol, weight_volume = 0.0;
    Teuchos::RCP<MeshIDs> ids = (*uspec)->second;

    for (MeshIDs::const_iterator c = ids->begin(); c != ids->end(); ++c) {
      if (*c < ncells_owned) {
        vol = mesh_->getCellVolume(*c);
        domain_volume_ += vol;
        weight_volume += vol * (*weight_)[0][*c];
      }
    }
    double result[2], tmp[2] = { domain_volume_, weight_volume };
    mesh_->getComm()->SumAll(tmp, result, 2);
    domain_volume_ = result[0];
    weight_volume = result[1];
    int nfun = (*uspec)->first->second->size();
    std::vector<double> val_vec(nfun);

    args[0] = t1;
    for (MeshIDs::const_iterator c = ids->begin(); c != ids->end(); ++c) {
      const AmanziGeometry::Point& xc = mesh_->getCellCentroid(*c);
      for (int i = 0; i != dim; ++i) args[i + 1] = xc[i];
      for (int i = 0; i < nfun; ++i) {
        val_vec[i] = (*(*uspec)->first->second)(args)[i] * (*weight_)[0][*c] / weight_volume;
      }
      value_[*c] = val_vec;
    }

    if (submodel_ == "integrated source") {
      args[0] = t0;
      for (MeshIDs::const_iterator c = ids->begin(); c != ids->end(); ++c) {
        const AmanziGeometry::Point& xc = mesh_->getCellCentroid(*c);
        for (int i = 0; i != dim; ++i) args[i + 1] = xc[i];
        for (int i = 0; i < nfun; ++i) {
          value_[*c][i] -= (*(*uspec)->first->second)(args)[i] * (*weight_)[0][*c] / weight_volume;
          value_[*c][i] *= dt;
        }
      }
    }
  }
}

} // namespace Amanzi

#endif
