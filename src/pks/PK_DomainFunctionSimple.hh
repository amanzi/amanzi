/*
  Process Kernels 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_PK_DOMAIN_FUNCTION_SIMPLE_HH_
#define AMANZI_PK_DOMAIN_FUNCTION_SIMPLE_HH_

#include <string>
#include <vector>

#include "Epetra_Vector.h"
#include "Teuchos_RCP.hpp"

#include "CommonDefs.hh"
#include "Mesh.hh"
#include "UniqueMeshFunction.hh"

namespace Amanzi {

template <class FunctionBase>
class PK_DomainFunctionSimple : public FunctionBase,
                                public Functions::UniqueMeshFunction {
 public:
  PK_DomainFunctionSimple(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) :
      UniqueMeshFunction(mesh) {};
  ~PK_DomainFunctionSimple() {};

  // member functions
  void Init(const Teuchos::ParameterList& plist);

  // required member functions
  virtual void Compute(double t0, double t1);
  virtual std::string name() { return "simple"; }

 protected:
  using FunctionBase::value_;

 private:
  std::string submodel_;
};


template <class FunctionBase>
void PK_DomainFunctionSimple<FunctionBase>::Init(const Teuchos::ParameterList& plist)
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
* Compute and distribute the result by Simple.
****************************************************************** */
template <class FunctionBase>
void PK_DomainFunctionSimple<FunctionBase>::Compute(double t0, double t1)
{
  if (unique_specs_.size() == 0) return;

  // create the input tuple (time + space)
  int dim = mesh_->space_dimension();
  std::vector<double> args(1 + dim);

  for (UniqueSpecList::const_iterator uspec = unique_specs_[AmanziMesh::CELL]->begin();
       uspec != unique_specs_[AmanziMesh::CELL]->end(); ++uspec) {

    args[0] = t1;
    Teuchos::RCP<MeshIDs> ids = (*uspec)->second;

    for (MeshIDs::const_iterator c = ids->begin(); c != ids->end(); ++c) {
      const AmanziGeometry::Point& xc = mesh_->cell_centroid(*c);
      for (int i = 0; i != dim; ++i) args[i + 1] = xc[i];

      // uspec->first is a RCP<Spec>, Spec's second is an RCP to the function.
      value_[*c] = (*(*uspec)->first->second)(args)[0];
    }
   
    if (submodel_ == "integrated source") {
      double dt = t1 - t0;
      if (dt > 0.0) dt = 1.0 / dt;
 
      args[0] = t0;
      for (MeshIDs::const_iterator c = ids->begin(); c != ids->end(); ++c) {
        const AmanziGeometry::Point& xc = mesh_->cell_centroid(*c);
        for (int i = 0; i != dim; ++i) args[i + 1] = xc[i];

        value_[*c] -= (*(*uspec)->first->second)(args)[0];
        value_[*c] *= dt;
       }
    }
  }
}

}  // namespace Amanzi

#endif
