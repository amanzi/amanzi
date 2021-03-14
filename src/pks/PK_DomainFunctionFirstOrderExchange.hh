/*
  Process Kernels 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (coonet@ornl.gov)
*/

#ifndef AMANZI_PK_DOMAIN_FUNCTION_FIRSTORDER_EXCHANGE_HH_
#define AMANZI_PK_DOMAIN_FUNCTION_FIRSTORDER_EXCHANGE_HH_

#include <string>
#include <vector>

#include "Epetra_Vector.h"
#include "Teuchos_RCP.hpp"

#include "CommonDefs.hh"
#include "Mesh.hh"
#include "UniqueMeshFunction.hh"
#include "DenseVector.hh"


namespace Amanzi {

template <class FunctionBase>
class PK_DomainFunctionFirstOrderExchange : public FunctionBase,
                                public Functions::UniqueMeshFunction {
 public:
  PK_DomainFunctionFirstOrderExchange(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                          const Teuchos::ParameterList& plist,
                          AmanziMesh::Entity_kind kind) :
    UniqueMeshFunction(mesh),
    FunctionBase(plist),
    kind_(kind) {
  };
  virtual ~PK_DomainFunctionFirstOrderExchange() = default;

  // member functions
  void Init(const Teuchos::ParameterList& plist, const std::string& keyword);

  virtual void set_state(const Teuchos::RCP<State>& S) { S_ = S; }

  // required member functions
  virtual void Compute(double t0, double t1);
  virtual std::string name() const { return "first order exchange"; }

 protected:
  using FunctionBase::value_;
  using FunctionBase::keyword_;

  Teuchos::RCP<const State> S_;

 private:
  AmanziMesh::Entity_kind kind_;
  std::string tcc_key_;
  std::string tcc_copy_;
};


/* ******************************************************************
* Initialization adds a single function to the list of unique specs.
****************************************************************** */
template <class FunctionBase>
void PK_DomainFunctionFirstOrderExchange<FunctionBase>::Init(
    const Teuchos::ParameterList& plist, const std::string& keyword)
{
  keyword_ = keyword;

  // get the model parameters
  Teuchos::ParameterList blist = plist.sublist("source function");
  tcc_key_ = Keys::readKey(blist, blist.get<std::string>("domain name", "domain"),
                           "total component concentration",
                           "total_component_concentration");
  tcc_copy_ = blist.get<std::string>("total component concentration copy", "default");
  
  // get and check the regions
  std::vector<std::string> regions =
      plist.get<Teuchos::Array<std::string> >("regions").toVector();

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
  Teuchos::RCP<Domain> domain = Teuchos::rcp(new Domain(regions, kind_));
  AddSpec(Teuchos::rcp(new Spec(domain, f)));
}


/* ******************************************************************
* Compute and distribute the result by FirstOrderExchange.
****************************************************************** */
template <class FunctionBase>
void PK_DomainFunctionFirstOrderExchange<FunctionBase>::Compute(double t0, double t1)
{
  if (unique_specs_.size() == 0) return;

  // create the input tuple (time + space)
  int dim = mesh_->space_dimension();
  std::vector<double> args(1 + dim);

  // get the tcc vector
  auto& tcc = *S_->GetFieldCopyData(tcc_key_, tcc_copy_)
              ->ViewComponent("cell", false);
  
  for (UniqueSpecList::const_iterator uspec = unique_specs_.at(kind_)->begin();
       uspec != unique_specs_.at(kind_)->end(); ++uspec) {

    args[0] = t1;
    Teuchos::RCP<MeshIDs> ids = (*uspec)->second;
    // uspec->first is a RCP<Spec>, Spec's second is an RCP to the function.
    int nfun = (*uspec)->first->second->size();
    std::vector<double> val_vec(nfun,0.);
    
    for (auto c = ids->begin(); c != ids->end(); ++c) {
      auto xc = PKUtils_EntityCoordinates(*c, kind_, *mesh_);

      for (int i = 0; i != dim; ++i) args[i + 1] = xc[i];
            
      // uspec->first is a RCP<Spec>, Spec's second is an RCP to the function.
      for (int i = 0; i < nfun; ++i) {
        val_vec[i] = -(*(*uspec)->first->second)(args)[i]*tcc[i][*c];
      }
      value_[*c] = val_vec;
    }
  }
}

}  // namespace Amanzi

#endif
