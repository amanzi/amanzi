/*
  Process Kernels 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (coonet@ornl.gov)
*/

#ifndef AMANZI_PK_DOMAIN_FUNCTION_SUBGRID_RETURN_HH_
#define AMANZI_PK_DOMAIN_FUNCTION_SUBGRID_RETURN_HH_

#include <string>
#include <vector>

#include "Epetra_Vector.h"
#include "Teuchos_RCP.hpp"

#include "CommonDefs.hh"
#include "DenseVector.hh"
#include "Mesh.hh"
#include "State.hh"
#include "Tag.hh"


namespace Amanzi {

template <class FunctionBase>
class PK_DomainFunctionSubgridReturn : public FunctionBase,
                                       public Functions::UniqueMeshFunction {
  
 public:
  PK_DomainFunctionSubgridReturn(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
          const Teuchos::ParameterList& plist) :
      FunctionBase(plist),
      UniqueMeshFunction(mesh) {};
  virtual ~PK_DomainFunctionSubgridReturn() = default;

  // member functions
  void Init(const Teuchos::ParameterList& plist, const std::string& keyword);

  // required member functions
  virtual void Compute(double t0, double t1);  
  virtual std::string name() const { return dset_; }
  virtual void set_state(const Teuchos::RCP<State>& S) { S_ = S; }

 protected:
  using FunctionBase::value_;
  using FunctionBase::keyword_;
  using Functions::UniqueMeshFunction::mesh_;
  Teuchos::RCP<const State> S_;

 private:
  Key field_out_suffix_, dset_; 
  Tag copy_field_out_key_; 

  Key saturation_key_, porosity_key_, molar_density_key_;
  Tag saturation_copy_, porosity_copy_, molar_density_copy_;
};


/* ******************************************************************
* Initialization
****************************************************************** */
template <class FunctionBase>
void PK_DomainFunctionSubgridReturn<FunctionBase>::Init(
    const Teuchos::ParameterList& plist, const std::string& keyword)
{
  keyword_ = keyword;

  // get and check the model parameters
  Teuchos::ParameterList blist = plist.sublist("source function");
  std::string domain = blist.get<std::string>("domain name", "surface"); 

  field_out_suffix_ = blist.get<std::string>("subgrid field suffix");
  copy_field_out_tag_ = make_tag(blist.get<std::string>("copy subgrid field", "default"));

  // check for prefixing -- mostly used in the multisubgrid models
  dset_ = blist.get<std::string>("subgrid domain set", "subgrid");

  // can "surface" prefix for this field be read automatically?? it is hardcoded now, or may be not?? 
  saturation_key_= Keys::readKey(blist, domain, "saturation liquid", "ponded_depth");
  saturation_copy_ = make_tag(blist.get<std::string>("saturation liquid copy", "default"));

  porosity_key_= Keys::readKey(blist, domain, "porosity", "porosity");
  porosity_copy_ = make_tag(blist.get<std::string>("porosity copy", "default"));
  
  molar_density_key_= Keys::readKey(blist, domain, "molar density", "molar_density_liquid");
  molar_density_copy_ = make_tag(blist.get<std::string>("molar density copy", "default"));

  // get and check the region
  auto regions = plist.get<Teuchos::Array<std::string> >("regions").toVector();

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
  Teuchos::RCP<Domain> domain = Teuchos::rcp(new Domain(regions, AmanziMesh::CELL));
  AddSpec(Teuchos::rcp(new Spec(domain, f)));
  
}


/* ******************************************************************
* Compute and distribute the result by SubgridReturn.
****************************************************************** */
template <class FunctionBase>
void PK_DomainFunctionSubgridReturn<FunctionBase>::Compute(double t0, double t1)
{
  if (unique_specs_.size() == 0) return;
  
  // create the input tuple (time + space)
  int dim = mesh_->space_dimension();
  std::vector<double> args(1 + dim);
  
  // get the map to convert to subgrid GID
  auto& map = mesh_->map(AmanziMesh::CELL, false);

  const auto& ws_ = *S_->Get<CompositeVector>(saturation_key_, saturation_copy_).ViewComponent("cell");
  const auto& phi_ = *S_->Get<CompositeVector>(porosity_key_, porosity_copy_).ViewComponent("cell");
  const auto& mol_dens_ = *S_->Get<CompositeVector>(molar_density_key_, molar_density_copy_).ViewComponent("cell");
  
  for (auto uspec = unique_specs_.at(AmanziMesh::CELL)->begin();
       uspec != unique_specs_.at(AmanziMesh::CELL)->end(); ++uspec) {
    
    args[0] = t1;
    Teuchos::RCP<MeshIDs> ids = (*uspec)->second;
    
    // uspec->first is a RCP<Spec>, Spec's second is an RCP to the function.
    int nfun = (*uspec)->first->second->size();
    
    // loop over all entities in the spec
    for (MeshIDs::const_iterator c = ids->begin(); c != ids->end(); ++c) {
      const AmanziGeometry::Point& xc = mesh_->cell_centroid(*c);
      
      // set the xyz coords
      for (int i = 0; i != dim; ++i) args[i + 1] = xc[i];
      
      // uspec->first is a RCP<Spec>, Spec's second is an RCP to the function.
      std::vector<double> alpha(nfun);
      
      for (int i = 0; i < nfun; ++i) {
        alpha[i] = (*(*uspec)->first->second)(args)[i];
      }
      
      // find the subgrid gid to be integrated
      auto gid = map.GID(*c);
      
      Key domain = Keys::getDomainInSet(name(), gid);
      
      Key var = Keys::getKey(domain,field_out_suffix_);
      
      // get the vector to be integrated
      const auto& vec_out = S_->Get<CompositeVector>(var, copy_field_out_tag_);
      const auto& vec_c = *vec_out.ViewComponent("cell", true);
      
      std::vector<double> val(nfun,0.);
      
      // DO THE INTEGRAL: currently omega_i = 1/cv_sg?
      int ncells_sg = vec_out.Mesh()->num_entities(AmanziMesh::CELL, AmanziMesh::Parallel_type::ALL);
      
      for (int c_sg = 0; c_sg! = ncells_sg; ++c_sg) {
        for (int k = 0; k != nfun; ++k) {
          val[k] += vec_c[k][c_sg] * alpha[k];
        }
      }
      
      for (int k=0; k!=nfun; ++k) {
        val[k] /= ncells_sg;
        val[k] *= ws_[k][*c] * phi_[k][*c] * mol_dens_[k][*c];
      }
      value_[*c] = val;
    }
  }
}
  
}  // namespace Amanzi

#endif
