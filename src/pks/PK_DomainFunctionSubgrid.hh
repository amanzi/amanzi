/*
  Process Kernels 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (coonet@ornl.gov)
*/

#ifndef AMANZI_PK_DOMAIN_FUNCTION_SUBGRID_HH_
#define AMANZI_PK_DOMAIN_FUNCTION_SUBGRID_HH_

#include <string>
#include <vector>

#include "Epetra_Vector.h"
#include "Teuchos_RCP.hpp"

#include "CommonDefs.hh"
#include "Mesh.hh"
#include "State.hh"
#include "DenseVector.hh"


namespace Amanzi {

template <class FunctionBase>
class PK_DomainFunctionSubgrid : public FunctionBase {
 public:
  PK_DomainFunctionSubgrid(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) :
      mesh_(mesh) {};
  virtual ~PK_DomainFunctionSubgrid() = default;

  // member functions
  void Init(const Teuchos::ParameterList& plist, const std::string& keyword,
            AmanziMesh::Entity_kind kind);

  // required member functions
  virtual void Compute(double t0, double t1);  
  virtual std::string name() const { return "hyporheic exchange"; }

  virtual void set_state(const Teuchos::RCP<State>& S) { S_ = S; }
  
 protected:
  using FunctionBase::value_;
  using FunctionBase::keyword_;

  Teuchos::RCP<const State> S_; 
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  
 private:
  AmanziMesh::Entity_ID entity_lid_;
  AmanziMesh::Entity_ID entity_gid_out_;

  std::string field_out_key_, copy_field_out_key_;
};


/* ******************************************************************
* Initialization adds a single function to the list of unique specs.
****************************************************************** */
template <class FunctionBase>
void PK_DomainFunctionSubgrid<FunctionBase>::Init(
    const Teuchos::ParameterList& plist, const std::string& keyword,
    AmanziMesh::Entity_kind region_kind)
{
  keyword_ = keyword;

  try {
    Teuchos::ParameterList blist = plist.sublist(keyword);
    field_out_key_ = blist.get<std::string>("field_out_key");
    if (blist.isParameter("copy_field_out_key"))
      copy_field_out_key_ = blist.get<std::string>("copy_field_out_key");
    else
      copy_field_out_key_ = "default";
  } catch (Errors::Message& msg) {
    Errors::Message m;
    m << "error in domain subgrid sublist : " << msg.what();
    Exceptions::amanzi_throw(m);
  }
    
  auto regions = plist.get<Teuchos::Array<std::string> >("regions");
  if (regions.size() != 1) {
    Errors::Message m;
    m << "Invalid \"regions\" list in processing subgrid source spec (must be of length 1).";
    Exceptions::amanzi_throw(m);
  }
  std::string region = regions[0];
  if (mesh_->valid_set_name(region, region_kind)) {
    AmanziMesh::Entity_ID_List id_list;
    mesh_->get_set_entities(region, region_kind, AmanziMesh::Parallel_type::ALL, &id_list);

    if (id_list.size() != 1) {
      Errors::Message m;
      m << "Invalid region in processing subgrid source spec (must be of length 1): \"" << region << "\"";
      Exceptions::amanzi_throw(m);
    }    
    entity_lid_ = id_list[0];

  } else {
    Errors::Message m;
    m << "Unknown region in processing subgrid source spec: \"" << region << "\"";
    Exceptions::amanzi_throw(m);
  }

  entity_gid_out_ = plist.get<int>("entity_gid_out");
}


/* ******************************************************************
* Compute and distribute the result by Subgrid.
****************************************************************** */
template <class FunctionBase>
void PK_DomainFunctionSubgrid<FunctionBase>::Compute(double t0, double t1)
{
  auto& vec_out = *S_->GetFieldCopyData(field_out_key_, copy_field_out_key_);
          
  const Epetra_MultiVector& field_out = *vec_out.ViewComponent("cell", true);
  int num_vec = field_out.NumVectors();
  std::vector<double> val(num_vec);

  AmanziMesh::Entity_ID entity_lid_out = vec_out.Mesh()->cell_map("false").LID(entity_gid_out_);
  AMANZI_ASSERT(entity_lid_out >= 0);
  for (int k=0; k<num_vec; ++k) val[k] = field_out[k][entity_lid_out];
  value_[entity_lid_] = val;
}

}  // namespace Amanzi

#endif
