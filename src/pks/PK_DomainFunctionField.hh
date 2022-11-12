/*
  Process Kernels 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef AMANZI_PK_DOMAIN_FUNCTION_FIELD_HH_
#define AMANZI_PK_DOMAIN_FUNCTION_FIELD_HH_

#include <string>
#include <vector>

#include "Epetra_Vector.h"
#include "Teuchos_RCP.hpp"

#include "CommonDefs.hh"
#include "Mesh.hh"
#include "UniqueMeshFunction.hh"
#include "Evaluator.hh"

#include "PK_Utils.hh"

namespace Amanzi {

template <class FunctionBase>
class PK_DomainFunctionField : public FunctionBase{
                             
 public:
  PK_DomainFunctionField(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                         const Teuchos::RCP<State>& S,
                          AmanziMesh::Entity_kind kind) :
    S_(S),
    mesh_(mesh),
    kind_(kind) {};

  PK_DomainFunctionField(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                         const Teuchos::RCP<State>& S,
                         const Teuchos::ParameterList& plist,
                         AmanziMesh::Entity_kind kind) :
    S_(S),
    mesh_(mesh),
    FunctionBase(plist),
    kind_(kind) {
  };

  ~PK_DomainFunctionField() {};

  typedef std::vector<std::string> RegionList;
  typedef std::pair<RegionList, AmanziMesh::Entity_kind> Domain;
  typedef std::set<AmanziMesh::Entity_ID> MeshIDs;


  // member functions
  void Init(const Teuchos::ParameterList& plist, const std::string& keyword);

  // required member functions
  virtual void Compute(double t0, double t1);
  virtual std::string name() const { return "field"; }

 protected:
  using FunctionBase::value_;
  using FunctionBase::keyword_;  
  Teuchos::RCP<State> S_;
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;

 private:
  std::string submodel_;
  Teuchos::RCP<MeshIDs> entity_ids_;
  AmanziMesh::Entity_kind kind_;
  Key field_key_;
  Tag tag_;
  Key component_key_;
};


/* ******************************************************************
* Initialization adds a single function to the list of unique specs.
****************************************************************** */
template <class FunctionBase>
void PK_DomainFunctionField<FunctionBase>::Init(
    const Teuchos::ParameterList& plist, const std::string& keyword)
{
  Errors::Message msg;

  keyword_ = keyword;

  submodel_ = "rate";
  if (plist.isParameter("submodel"))
    submodel_ = plist.get<std::string>("submodel");
  std::vector<std::string> regions = plist.get<Teuchos::Array<std::string> >("regions").toVector();
  
  // Teuchos::RCP<Amanzi::MultiFunction> f;
  try {
    Teuchos::ParameterList flist = plist.sublist(keyword);
    field_key_ = flist.get<std::string>("field key");
    tag_ = Keys::readTag(flist, "tag");
    component_key_ = flist.get<std::string>("component", "cell");;    
  } catch (Errors::Message& msg) {
    Errors::Message m;
    m << "error in source sublist : " << msg.what();
    Exceptions::amanzi_throw(m);
  }

  // Add this source specification to the domain function.
  Teuchos::RCP<Domain> domain = Teuchos::rcp(new Domain(regions, kind_));
  //AddSpec(Teuchos::rcp(new Spec(domain, f)));

  entity_ids_ = Teuchos::rcp(new MeshIDs());
  AmanziMesh::Entity_kind kind = domain->second;

  for (auto region = domain->first.begin(); region != domain->first.end(); ++region) {
    if (mesh_->valid_set_name(*region, kind)) {
      AmanziMesh::Entity_ID_List id_list;
      mesh_->get_set_entities(*region, kind, AmanziMesh::Parallel_type::OWNED, &id_list);
      entity_ids_->insert(id_list.begin(), id_list.end());
    } else {
      msg << "Unknown region in processing coupling source: name=" << *region
          << ", kind=" << kind << "\n";
      Exceptions::amanzi_throw(msg);
    }
  }


  
}


/* ******************************************************************
* Compute and distribute the result by Simple.
****************************************************************** */
template <class FunctionBase>
void PK_DomainFunctionField<FunctionBase>::Compute(double t0, double t1)
{
  // create the input tuple (time + space)
  int dim = mesh_->space_dimension();

  if (S_->HasEvaluator(field_key_, tag_)){
    S_->GetEvaluator(field_key_, tag_).Update(*S_, field_key_);    
  }  
  
  const auto& field_vec = *S_->Get<CompositeVector>(field_key_, tag_).ViewComponent("cell", false);
  
  int nvalues = field_vec.NumVectors();  
  std::vector<double> val_vec(nvalues);  

  for (auto c : *entity_ids_) {
      for (int i = 0; i < nvalues; ++i) {
        val_vec[i] = field_vec[i][c];
      }
      value_[c] = val_vec;
  }
}

}  // namespace Amanzi

#endif
