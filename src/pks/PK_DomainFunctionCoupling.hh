/*
  Process Kernels

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Author: Daniil Svyatsky (dasvyat@lanl.gov)

  This provide coupling of fields located on matching manifold and
  space meshes. For the space, the coupling creates a list of boundary
  conditions. For the manifold, the coupling creates a list of sources.

  Typically, the manifold mesh provides map manifold (cell) -> space (face).
  In space, we need the reverse map.

  The following parameter names were changed:
    "flux_key" -> "flux key"

    "field_in_key" -> 
    "copy_field_in_key" -> 

    "field_out_key" -> "field key"
    "copy_field_out_key" -> "field copy key"

  There are three submodels.
  (A) submodel="rate". The computed data are given by formulas: 

      value[i][c] += flux[f] * field[i][cc] / |c|
      value[N][c] -= dt * flux[f] 

  where cc is the space cell attached to face f, and N is the auxiliary 
  value added to the result. Note that internal face f (resp., boundary 
  face f) is shared by two (resp. one) space cells.

  (B) submodel="field". The computed data are given by the formula:

      value[i][f] = F(field[i][c]); 

  where F is function.

  (C) submodel="conserved quantity". Not used in Amanzi.
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
#include "FunctionFactory.hh"
#include "Mesh.hh"
#include "State.hh"
#include "UniqueLocalIndex.hh"

namespace Amanzi {

template <class FunctionBase>
class PK_DomainFunctionCoupling : public FunctionBase {
 public:
  PK_DomainFunctionCoupling(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                            const Teuchos::RCP<const State>& S) :
      mesh_(mesh), S_(S) {};
  virtual ~PK_DomainFunctionCoupling() = default;

  typedef std::vector<std::string> RegionList;
  typedef std::pair<RegionList, AmanziMesh::Entity_kind> Domain;
  typedef std::set<AmanziMesh::Entity_ID> MeshIDs;

  // member functions
  void Init(const Teuchos::ParameterList& plist, const std::string& keyword,
            AmanziMesh::Entity_kind kind);

  // required member functions
  virtual void Compute(double t0, double t1) override;
  virtual std::string name() const override { return "domain coupling"; }

 protected:
  using FunctionBase::value_;
  using FunctionBase::linear_term_;
  using FunctionBase::keyword_;

  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  Teuchos::RCP<const State> S_;

 private:
  std::string submodel_;

  std::string flux_key_, copy_flux_key_;
  std::string field_cons_key_, copy_field_cons_key_;
  std::string field_out_key_, copy_field_out_key_;
  std::string field_in_key_, copy_field_in_key_;

  std::unique_ptr<Function> f_;

  Teuchos::RCP<MeshIDs> entity_ids_;
  std::map<AmanziMesh::Entity_ID, AmanziMesh::Entity_ID> reverse_map_;
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

  Teuchos::ParameterList slist;
  try {
    slist = plist.sublist(keyword_);
  } catch (Errors::Message& msg) {
    Errors::Message m;
    m << "error in domain coupling sublist : " << msg.what();
    Exceptions::amanzi_throw(m);
  }

  // deprecated
  /*
  std::vector<std::string> deprecated = { "flux_key", "copy_flux_key",
                                          "field_in_key", "copy_field_in_key",
                                          "field_out_key", "copy_field_out_key" };

  for (auto name : deprecated) {
    if (slist.isParameter(name) {
      Errors::Message msg;
      msg << "deprecated name : " << name << ", see PK_DomainFunctionCoupling.hh";
      Exceptions::amanzi_throw(msg);
    }
  }
  */

  // get keys of owned (in) and exterior (out) fields
  if (submodel_ == "rate") {
    field_in_key_ = slist.get<std::string>("field_in_key");
    copy_field_in_key_ = slist.get<std::string>("copy_field_in_key", "default");

    flux_key_ = slist.get<std::string>("flux_key");
    copy_flux_key_ = slist.get<std::string>("copy_flux_key", "default");

  } else if (submodel_ == "rate_Darcy") {
    field_in_key_ = slist.get<std::string>("field_in_key");
    copy_field_in_key_ = slist.get<std::string>("copy_field_in_key", "default");

    flux_key_ = slist.get<std::string>("flux_key");
    copy_flux_key_ = slist.get<std::string>("copy_flux_key", "default");

  } else if (submodel_ == "field") {
    field_in_key_ = slist.get<std::string>("field_in_key");
    copy_field_in_key_ = slist.get<std::string>("copy_field_in_key", "default");

  } else if (submodel_ == "conserved quantity") {
    field_cons_key_ = slist.get<std::string>("conserved_quantity_key");
    copy_field_cons_key_ = slist.get<std::string>("copy_conserved_quantity_key", "default");

  } else {
    Errors::Message m;
    m << "unknown DomainFunctionCoupling submodel \"" << submodel_ << "\", valid are: \"field\", \"rate\", and \"conserved quantity\"";
    Exceptions::amanzi_throw(m);
  }

  field_out_key_ = slist.get<std::string>("field_out_key");
  copy_field_out_key_ = slist.get<std::string>("copy_field_out_key", "default");

  // create a list of domain ids
  RegionList regions = plist.get<Teuchos::Array<std::string> >("regions").toVector();
  Teuchos::RCP<Domain> domain= Teuchos::rcp(new Domain(regions, region_kind));

  entity_ids_ = Teuchos::rcp(new MeshIDs());
  AmanziMesh::Entity_kind kind = domain->second;

  for (auto region = domain->first.begin(); region != domain->first.end(); ++region) {
    if (mesh_->valid_set_name(*region, kind)) {
      AmanziMesh::Entity_ID_List id_list;
      mesh_->get_set_entities(*region, kind, AmanziMesh::Parallel_type::OWNED, &id_list);
      entity_ids_->insert(id_list.begin(), id_list.end());
    } else {
      std::stringstream m;
      m << "Unknown region in processing coupling source: name=" << *region
        << ", kind=" << kind << "\n";
      Errors::Message message(m.str());
      Exceptions::amanzi_throw(message);
    }
  }

  // create a function
  if (slist.isSublist("function")) {
    auto flist = slist.sublist("function");
    FunctionFactory factory;
    f_ = factory.Create(flist);
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
    const auto& flux = *S_->GetFieldCopyData(flux_key_, copy_flux_key_)->ViewComponent("face", true);
    const auto& field_out = *S_->GetFieldCopyData(field_out_key_, copy_field_out_key_)->ViewComponent("cell", true);

    const auto& flux_map = S_->GetFieldData(flux_key_)->Map().Map("face", true);
    const auto& source_mesh = S_->GetMesh(Keys::getDomain(field_in_key_));

    S_->GetFieldCopyData(field_out_key_, copy_field_out_key_)->ScatterMasterToGhosted("cell");

    int num_vec = field_out.NumVectors();

    AmanziMesh::Entity_ID_List cells, faces;
    std::vector<int> dirs;
    auto mesh_out = S_->GetFieldData(field_out_key_)->Mesh();

    // loop over cells on the manifold
    for (auto c : *entity_ids_) {
      AmanziMesh::Entity_ID f = mesh_->entity_get_parent(AmanziMesh::CELL, c);

      mesh_out->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);

      if (cells.size() != flux_map->ElementSize(f)) {
        std::stringstream m;
        m << "Number of flux DOFs doesn't equal to the number of cells sharing a face: "
          << cells.size() << " != " << flux_map->ElementSize(f) << std::endl;
        Errors::Message message(m.str());
        Exceptions::amanzi_throw(message);
      }

      double linear(0.0);
      std::vector<double> val(num_vec + 1, 0.0);  // extra is for water
      int pos = Operators::UniqueIndexFaceToCells(*mesh_out, f, cells[0]);

      for (int j = 0; j != cells.size(); ++j) {
        mesh_out->cell_get_faces_and_dirs(cells[j], &faces, &dirs);

        for (int i = 0; i < faces.size(); i++) {
          if (f == faces[i]) {
            int g = flux_map->FirstPointInElement(f);
            double fln = flux[0][g + (pos + j)%2] * dirs[i];

std::cout << f << " " << mesh_out->face_centroid(f) << fln << std::endl;
            if (fln >= 0) {
              linear += fln;
              for (int k = 0; k < num_vec; ++k) {
                // flux is water flux * Conc
                val[k] += field_out[k][cells[j]] * fln / source_mesh->cell_volume(c);
              }
            } else {
              val[num_vec] += -fln * (t1 - t0);
            }
            break;
          }
        }
      }
      value_[c] = val;
      linear_term_[c] = linear;
    }

  } else if (submodel_ == "rate_Darcy") {
    const auto& flux = *S_->GetFieldCopyData(flux_key_, copy_flux_key_)->ViewComponent("face", true);
    const auto& field_out = *S_->GetFieldCopyData(field_out_key_, copy_field_out_key_)->ViewComponent("cell", true);
    const auto& field_in = *S_->GetFieldCopyData(field_in_key_, copy_field_in_key_)->ViewComponent("cell", true);

    const auto& flux_map = S_->GetFieldData(flux_key_)->Map().Map("face", true);
    const auto& source_mesh = S_->GetMesh(Keys::getDomain(field_in_key_));

    S_->GetFieldCopyData(field_out_key_, copy_field_out_key_)->ScatterMasterToGhosted("cell");
    S_->GetFieldCopyData(field_in_key_, copy_field_in_key_)->ScatterMasterToGhosted("cell");

    if (field_in.NumVectors() != field_out.NumVectors()) {
      std::stringstream m;
      m << "Mismatch of numbers of vectors in coupled fields.\n";
      Errors::Message message(m.str());
      Exceptions::amanzi_throw(message);
    }
    int num_vec = field_in.NumVectors();

    AmanziMesh::Entity_ID_List cells, faces;
    std::vector<int> dirs;
    auto mesh_out = S_->GetFieldData(field_out_key_)->Mesh();

    // loop over cells on the manifold
    for (auto c : *entity_ids_) {
      AmanziMesh::Entity_ID f = mesh_->entity_get_parent(AmanziMesh::CELL, c);

      mesh_out->face_get_cells(f, AmanziMesh::Parallel_type::ALL, &cells);

      if (cells.size() != flux_map->ElementSize(f)) {
        std::stringstream m;
        m << "Number of flux DOFs doesn't equal to the number of cells sharing a face: "
          << cells.size() << " != " << flux_map->ElementSize(f) << std::endl;
        Errors::Message message(m.str());
        Exceptions::amanzi_throw(message);
      }

      double linear(0.0);
      std::vector<double> val(num_vec + 1, 0.0);  // val(end) is for water
      int pos = Operators::UniqueIndexFaceToCells(*mesh_out, f, cells[0]);

      for (int j = 0; j != cells.size(); ++j) {
        mesh_out->cell_get_faces_and_dirs(cells[j], &faces, &dirs);

        for (int i = 0; i < faces.size(); i++) {
          if (f == faces[i]) {
            int g = flux_map->FirstPointInElement(f);
            double fln = flux[0][g + (pos + j)%2] * dirs[i];
            
            linear += fln;
            for (int k = 0; k < num_vec; ++k) {
              // flux is water flux * Conc
              val[k] += field_out[k][cells[j]] * fln / source_mesh->cell_volume(c);
            }
              
            val[num_vec] += fln * (t1 - t0);
            
            break;
          }
        }
      }
      value_[c] = val;
      linear_term_[c] = linear;
    }

  } else if (submodel_ == "field") {
    // create reserse map from space (face) onto manifold (cell)
    auto mesh_out = S_->GetFieldData(field_out_key_)->Mesh();

    if (mesh_->space_dimension() == mesh_->manifold_dimension() && reverse_map_.size() == 0) {
      const Epetra_Map& cell_map = mesh_out->cell_map(true);
      for (int c = 0 ; c < cell_map.NumMyElements(); ++c) {
        AmanziMesh::Entity_ID f = mesh_out->entity_get_parent(AmanziMesh::CELL, c);
        reverse_map_[f] = c;
      }
    }

    const auto& field_out = *S_->GetFieldCopyData(field_out_key_, copy_field_out_key_)->ViewComponent("cell", true);

    int num_vec = field_out.NumVectors();
    std::vector<double> val(num_vec);

    // Loop over faces (owned + ghosted) in the space restricted to the manifold.
    // The set of these faces could be bigger then the set of manifold cells (owned + ghosted)
    std::vector<double> args(1);

    for (auto f = entity_ids_->begin(); f != entity_ids_->end(); ++f) {
      auto it = reverse_map_.find(*f);
      if (it == reverse_map_.end()) continue;

      int c = it->second;
      for (int k = 0; k < num_vec; ++k) val[k] = field_out[k][c]; 

      if (f_ != nullptr) {
        for (int k = 0; k < num_vec; ++k) {
          args[0] = val[k];
          val[k] = (*f_)(args);
        }
      }
      value_[*f] = val;     
    }

  } else if (submodel_ == "conserved quantity") {
    // create reverse map from space (face) onto manifold (cell)
    auto mesh_out = S_->GetFieldData(field_out_key_)->Mesh();

    if (reverse_map_.size() == 0 && mesh_->space_dimension() == mesh_->manifold_dimension()) {
      const Epetra_Map& cell_map = mesh_out->cell_map(true);
      for (int c = 0 ; c < cell_map.NumMyElements(); ++c) {
        AmanziMesh::Entity_ID f = mesh_out->entity_get_parent(AmanziMesh::CELL, c);
        reverse_map_[f] = c;
      }
    }

    const auto& field_out = *S_->GetFieldCopyData(field_out_key_, copy_field_out_key_)->ViewComponent("cell", true);
    const auto& field_cons = *S_->GetFieldCopyData(field_cons_key_, copy_field_cons_key_)->ViewComponent("cell", true);

    int num_vec = field_out.NumVectors();
    std::vector<double> val(num_vec);
    AMANZI_ASSERT(num_vec + 2 == field_cons.NumVectors());

    // Loop over faces (owned + ghosted) in the space restricted to the manifold.
    // The set of these faces could be bigger then the set of manifold cells (owned + ghosted)
    for (auto f : *entity_ids_) {
      auto it = reverse_map_.find(f);
      if (it == reverse_map_.end()) continue;
      int sc = it->second;

      // accept it all
      AmanziMesh::Entity_ID_List cells, faces;
      std::vector<int> dirs;
      mesh_->face_get_cells(f, AmanziMesh::Parallel_type::OWNED, &cells);
      AMANZI_ASSERT(cells.size() == 1);

      mesh_->cell_get_faces_and_dirs(cells[0], &faces, &dirs);
      auto i = std::find(faces.begin(), faces.end(), f) - faces.begin();
      AMANZI_ASSERT(i < dirs.size());

      for (int k = 0; k < num_vec; ++k) {
        if (field_cons[num_vec][sc] > 0.) {
          val[k] = field_cons[k][sc] / field_cons[num_vec][sc];
        } else {
          val[k] = field_out[k][sc];
        }
      }
      value_[f] = val;
    }
  }
}

}  // namespace Amanzi

#endif
