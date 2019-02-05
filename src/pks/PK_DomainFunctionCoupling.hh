/*
  Process Kernels 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Daniil Svyatsky (dasvyat@lanl.gov)
*/

#ifndef AMANZI_PK_DOMAIN_FUNCTION_COUPLING_HH_
#define AMANZI_PK_DOMAIN_FUNCTION_COUPLING_HH_

#include <string>
#include <vector>

// TPLs
#include "Epetra_Vector.h"
#include "Teuchos_RCP.hpp"

// Amanzi
#include "CommonDefs.hh"
#include "DenseVector.hh"
#include "Mesh.hh"
#include "State.hh"

namespace Amanzi {

template <class FunctionBase>
class PK_DomainFunctionCoupling : public FunctionBase {
 public:
  PK_DomainFunctionCoupling(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh) :
      mesh_(mesh) {};
  virtual ~PK_DomainFunctionCoupling() = default;

  typedef std::vector<std::string> RegionList;
  typedef std::pair<RegionList, AmanziMesh::Entity_kind> Domain;
  typedef std::set<AmanziMesh::Entity_ID> MeshIDs;

  // member functions
  void Init(const Teuchos::ParameterList& plist, const std::string& keyword,
            AmanziMesh::Entity_kind kind);

  // required member functions
  virtual void Compute(double t0, double t1);
  virtual std::string name() const { return "domain coupling"; }
  virtual void set_state(const Teuchos::RCP<State>& S) { S_ = S; }

 protected:
  using FunctionBase::value_;
  using FunctionBase::keyword_;

  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  Teuchos::RCP<const State> S_;
  
 private:
  std::string submodel_;
  std::string flux_key_, copy_flux_key_;
  std::string field_out_key_, copy_field_out_key_;
  std::string field_in_key_, copy_field_in_key_;
  Teuchos::RCP<MeshIDs> entity_ids_; 
};


/* ******************************************************************
* Initialization adds a single function to the list of unique specs.
****************************************************************** */
template <class FunctionBase>
void PK_DomainFunctionCoupling<FunctionBase>::Init(
    const Teuchos::ParameterList& plist, const std::string& keyword,
    AmanziMesh::Entity_kind region_kind)
{
  keyword_ = keyword;

  if (plist.isParameter("submodel")) {
    submodel_ = plist.get<std::string>("submodel");
  } else {
    Errors::Message m;
    m << "Domain coupling requires parameter \"submodel\"";
    Exceptions::amanzi_throw(m);
  }

  if (submodel_ == "rate") {
    try {
      Teuchos::ParameterList slist = plist.sublist("sink");
      if (slist.isParameter("field_in_key"))
        field_in_key_ = slist.get<std::string>("field_in_key");
      copy_field_in_key_ = slist.get<std::string>("copy_field_in_key", "default");
      if (slist.isParameter("field_out_key"))
        field_out_key_ = slist.get<std::string>("field_out_key");
      copy_field_out_key_ = slist.get<std::string>("copy_field_out_key", "default");
      if (slist.isParameter("flux_key"))
        flux_key_ = slist.get<std::string>("flux_key");
      copy_flux_key_ = slist.get<std::string>("copy_flux_key", "default");
    } catch (Errors::Message& msg) {
      Errors::Message m;
      m << "error in domain coupling sublist : " << msg.what();
      Exceptions::amanzi_throw(m);
    }
  } else if (submodel_ == "field") {
    try {
      Teuchos::ParameterList blist = plist.sublist("boundary concentration");
      if (blist.isParameter("field_out_key")) {
        field_out_key_ = blist.get<std::string>("field_out_key");
        copy_field_out_key_ = blist.get<std::string>("copy_field_out_key", "default");
      }
    } catch (Errors::Message& msg) {
      Errors::Message m;
      m << "error in domain coupling sublist : " << msg.what();
      Exceptions::amanzi_throw(m);
    }
  }
    
  RegionList regions = plist.get<Teuchos::Array<std::string> >("regions").toVector();
  // Add this source specification to the domain function.
  Teuchos::RCP<Domain> domain= Teuchos::rcp(new Domain(regions, region_kind));

  entity_ids_ = Teuchos::rcp(new MeshIDs());
  AmanziMesh::Entity_kind kind = domain->second;
  for (RegionList::const_iterator region = domain->first.begin();
      region != domain->first.end(); ++region) {
    // Get the ids from the mesh by region name and entity kind.
    if (mesh_->valid_set_name(*region, kind)) {
      AmanziMesh::Entity_ID_List id_list;
      mesh_->get_set_entities(*region, kind, AmanziMesh::Parallel_type::ALL, &id_list);
      entity_ids_->insert(id_list.begin(), id_list.end());
    } else {
      std::stringstream m;
      m << "Unknown region in processing coupling source spec: \"" << *region << "\"";
      Errors::Message message(m.str());
      Exceptions::amanzi_throw(message);
    }
  }
}


/* ******************************************************************
* Compute and distribute the result by Coupling.
****************************************************************** */
template <class FunctionBase>
void PK_DomainFunctionCoupling<FunctionBase>::Compute(double t0, double t1)
{
  // create the input tuple (time + space)
  if (submodel_ == "rate") {
    const Epetra_MultiVector& flux = 
        *S_->GetFieldCopyData(flux_key_, copy_flux_key_)->ViewComponent("face", true);

    Teuchos::RCP<const Epetra_BlockMap> flux_map = S_->GetFieldData(flux_key_)->Map().Map("face", true);

    S_->GetFieldCopyData(field_out_key_, copy_field_out_key_)->ScatterMasterToGhosted("cell");    
    const Epetra_MultiVector& field_out = 
        *S_->GetFieldCopyData(field_out_key_, copy_field_out_key_)->ViewComponent("cell", true);


    S_->GetFieldCopyData(field_in_key_, copy_field_in_key_)->ScatterMasterToGhosted("cell"); 
    const Epetra_MultiVector& field_in = 
        *S_->GetFieldCopyData(field_in_key_, copy_field_in_key_)->ViewComponent("cell", true);

    if (field_in.NumVectors() != field_out.NumVectors()) {
      std::stringstream m;
      m << "Mismatch of vector sizes in domain couplong function.\n";
      Errors::Message message(m.str());
      Exceptions::amanzi_throw(message);
    }
    int num_vec = field_in.NumVectors();


    Teuchos::RCP<const AmanziMesh::Mesh> mesh_out = S_->GetFieldData(field_out_key_)->Mesh();
    AmanziMesh::Entity_ID_List cells, faces;
    std::vector<int> dirs;
    const Epetra_Map& cell_map = mesh_out -> cell_map(true);

    
    for (MeshIDs::const_iterator c = entity_ids_->begin(); c != entity_ids_->end(); ++c) {

      AmanziMesh::Entity_ID f = mesh_->entity_get_parent(AmanziMesh::CELL, *c);

      mesh_out->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);

      // if (f==202)
      //   std::cout<<"Enter\n";
    
      if (cells.size() != flux_map -> ElementSize(f)) {
        std::stringstream m;
        m << "Number of flux DOF doesn't equal to the number of cell sharing the interface: cell_ids\"" << cells[0] << "\"";
        Errors::Message message(m.str());
        Exceptions::amanzi_throw(message);
      }

      std::vector<double> val(num_vec, 0);
      int pos = 0;
      if (cells.size() == 2) {
        pos = (cell_map.GID(cells[0]) < cell_map.GID(cells[1])) ? 0 : 1;
      }
      
      for (int j=0; j!=cells.size(); ++j) {

        mesh_out->cell_get_faces_and_dirs(cells[j], &faces, &dirs);

        for (int i = 0; i < faces.size(); i++) {
          if (f == faces[i]) {
            int f_loc_id = flux_map -> FirstPointInElement(f);            
            double fln = flux[0][f_loc_id + (pos + j)%2]*dirs[i];
            
            if (fln >= 0) {        
              for (int k=0; k<num_vec; ++k) {
                val[k] += field_out[k][cells[j]] * fln;
              }
            } else if (fln < 0) {       
              for (int k=0; k<num_vec; ++k) {
                val[k] += field_in[k][*c] * fln;
              }
            }
            break;
          }
        }
      }
      value_[*c] = val; 
    }
  } else if (submodel_ == "field") {
    const Epetra_MultiVector& field_out = 
        *S_->GetFieldCopyData(field_out_key_, copy_field_out_key_)->ViewComponent("cell", true);

    int num_vec = field_out.NumVectors();
    std::vector<double> val(num_vec);

    ///////CHECK!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!????????????
    int i(0);
    for (MeshIDs::const_iterator c = entity_ids_->begin(); c != entity_ids_->end(); ++c) {
      for (int k=0; k<num_vec; ++k) val[k] = field_out[k][i];            
      value_[*c] = val;     
      i++;
    }
  }
}

}  // namespace Amanzi

#endif
