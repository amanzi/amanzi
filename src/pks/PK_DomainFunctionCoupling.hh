/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatsky (dasvyat@lanl.gov)
*/

/*!

This provide coupling of fields located on matching manifold and
space meshes. For the space, the coupling creates a list of boundary
conditions. For the manifold, the coupling creates a list of sources.

Typically, the manifold mesh provides map manifold (cell) -> space (face).
In space, we need the reverse map.

The following parameter names were changed:
    "flux_key" -> "flux key"
    "copy_flux_tag" -> "flux copy key"

    "field_in_key" -> removed
    "copy_field_in_key" -> removed

    "field_out_key" -> "external field key"
    "copy_field_out_tag" -> "external field copy key"

    "conserved_quantity_key" -> "conserved quantity key"
    "copy_conserved_quantity_key" -> "conserved quantity copy key"

There are three submodels.
(A) submodel="rate". The computed data are given by formulas:

      value[i][c] += flux[f] * external_field[i][cc] / V_c

      value[N][c] -= dt * flux[f]

where cc is the space cell attached to face f, and N is the auxiliary
value added to the result. Note that internal face f (resp., boundary
face f) is shared by two (resp. one) space cells.

  (B) submodel="field". The computed data are given by the formula:

      value[i][f] = external_field[i][c];

  (C) submodel="conserved quantity". Not used in Amanzi.

  (D) submodel="".

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
#include "UniqueLocalIndex.hh"

namespace Amanzi {

template <class FunctionBase>
class PK_DomainFunctionCoupling : public FunctionBase {
 public:
  PK_DomainFunctionCoupling(const Teuchos::RCP<const AmanziMesh::Mesh>& mesh,
                            const Teuchos::RCP<const State>& S)
    : mesh_(mesh), S_(S){};
  virtual ~PK_DomainFunctionCoupling() = default;

  typedef std::vector<std::string> RegionList;
  typedef std::pair<RegionList, AmanziMesh::Entity_kind> Domain;
  typedef std::set<AmanziMesh::Entity_ID> MeshIDs;

  // member functions
  void Init(const Teuchos::ParameterList& plist,
            const std::string& keyword,
            AmanziMesh::Entity_kind kind);

  // required member functions
  virtual void Compute(double t0, double t1) override;
  virtual std::string name() const override { return "domain coupling"; }

 protected:
  using FunctionBase::value_;
  using FunctionBase::linear_term_;
  using FunctionBase::keyword_;

  Teuchos::RCP<const AmanziMesh::Mesh> mesh_, mesh_out_;
  Teuchos::RCP<const State> S_;

 private:
  std::string submodel_;

  Key flux_key_;
  Key field_cons_key_, field_out_key_;
  Tag copy_flux_tag_, copy_field_out_tag_, copy_field_cons_tag_;

  Teuchos::RCP<MeshIDs> entity_ids_;
  std::map<AmanziMesh::Entity_ID, AmanziMesh::Entity_ID> reverse_map_;
};


/* ******************************************************************
* Initialization adds a single function to the list of unique specs.
****************************************************************** */
template <class FunctionBase>
void
PK_DomainFunctionCoupling<FunctionBase>::Init(const Teuchos::ParameterList& plist,
                                              const std::string& keyword,
                                              AmanziMesh::Entity_kind region_kind)
{
  Errors::Message msg;

  keyword_ = keyword;

  if (plist.isParameter("submodel")) {
    submodel_ = plist.get<std::string>("submodel");
  } else {
    msg << "Domain coupling requires parameter \"submodel\"";
    Exceptions::amanzi_throw(msg);
  }

  Teuchos::ParameterList slist;
  try {
    slist = plist.sublist(keyword_);
  } catch (Errors::Message& msg2) {
    msg << "error in domain coupling sublist : " << msg2.what();
    Exceptions::amanzi_throw(msg);
  }

  //  deprecated
  std::vector<std::string> deprecated = {
    "flux_key",      "copy_flux_tag",      "field_in_key",           "copy_field_in_key",
    "field_out_key", "copy_field_out_tag", "conserved_quantity_key", "copy_conserved_quantity_key"
  };
  std::vector<std::string> newstyle = { "flux key",
                                        "flux copy key",
                                        "removed1",
                                        "removed2",
                                        "external field key",
                                        "external field copy key",
                                        "conserved quantity key",
                                        "conserved quantity copy key" };

  for (int i = 0; i < deprecated.size(); ++i) {
    std::string name = deprecated[i];
    if (slist.isParameter(name)) {
      slist.set<std::string>(newstyle[i], slist.get<std::string>(name));
      /*
      Errors::Message msg;
      msg << "deprecated name : " << name << ", see PK_DomainFunctionCoupling.hh";
      Exceptions::amanzi_throw(msg);
      */
    }
  }

  // get keys of owned (in) and exterior (out) fields
  if (submodel_ == "rate" || submodel_ == "flux exchange") {
    flux_key_ = slist.get<std::string>("flux key");
    copy_flux_tag_ = Keys::readTag(slist, "flux copy key");

  } else if (submodel_ == "conserved quantity") {
    field_cons_key_ = slist.get<std::string>("conserved quantity key");
    copy_field_cons_tag_ = Keys::readTag(slist, "conserved quantity copy key");

  } else if (submodel_ == "field") {
    // pass

  } else {
    msg << "unknown DomainFunctionCoupling submodel \"" << submodel_
        << "\", valid option: \"field\", \"rate\", \"flux exchange\", and \"conserved quantity\"";
    Exceptions::amanzi_throw(msg);
  }

  field_out_key_ = slist.get<std::string>("external field key");
  copy_field_out_tag_ = Keys::readTag(slist, "external field copy key");

  // create a list of domain ids
  RegionList regions = plist.get<Teuchos::Array<std::string>>("regions").toVector();
  Teuchos::RCP<Domain> domain = Teuchos::rcp(new Domain(regions, region_kind));

  entity_ids_ = Teuchos::rcp(new MeshIDs());
  AmanziMesh::Entity_kind kind = domain->second;

  for (auto region = domain->first.begin(); region != domain->first.end(); ++region) {
    if (mesh_->isValidSetName(*region, kind)) {
      auto id_list = mesh_->getSetEntities(*region, kind, AmanziMesh::Parallel_type::OWNED);
      entity_ids_->insert(id_list.begin(), id_list.end());
    } else {
      msg << "Unknown region in processing coupling source: name=" << *region << ", kind=" << kind
          << "\n";
      Exceptions::amanzi_throw(msg);
    }
  }
}


/* ******************************************************************
* Compute and distribute the result by Coupling.
****************************************************************** */
template <class FunctionBase>
void
PK_DomainFunctionCoupling<FunctionBase>::Compute(double t0, double t1)
{
  Errors::Message msg;

  // compute reverse map: should we move this to Init() ???
  mesh_out_ = S_->Get<CompositeVector>(field_out_key_, copy_field_out_tag_).Mesh();

  if (submodel_ == "field" || submodel_ == "conserved quantity") {
    if (mesh_->getSpaceDimension() == mesh_->getManifoldDimension() && reverse_map_.size() == 0) {
      const Epetra_Map& cell_map = mesh_out_->getMap(AmanziMesh::Entity_kind::CELL,true);
      for (int c = 0; c < cell_map.NumMyElements(); ++c) {
        AmanziMesh::Entity_ID f = mesh_out_->getEntityParent(AmanziMesh::Entity_kind::CELL, c);
        reverse_map_[f] = c;
      }
    }
  }

  // create the input tuple (time + space)
  if (submodel_ == "rate") {
    const auto& flux =
      *S_->Get<CompositeVector>(flux_key_, copy_flux_tag_).ViewComponent("face", true);
    const auto& flux_map =
      S_->Get<CompositeVector>(flux_key_, copy_flux_tag_).Map().Map("face", true);

    const auto& field_out =
      *S_->Get<CompositeVector>(field_out_key_, copy_field_out_tag_).ViewComponent("cell", true);
    S_->Get<CompositeVector>(field_out_key_, copy_field_out_tag_).ScatterMasterToGhosted("cell");

    int num_vec = field_out.NumVectors();

    AmanziMesh::Entity_ID_List cells;

    // loop over cells on the manifold
    for (auto c : *entity_ids_) {
      AmanziMesh::Entity_ID f = mesh_->getEntityParent(AmanziMesh::Entity_kind::CELL, c);

      cells = mesh_out_->getFaceCells(f, AmanziMesh::Parallel_type::ALL);

      if (cells.size() != flux_map->ElementSize(f)) {
        msg << "Number of flux DOFs doesn't equal to the number of cells sharing a face: "
            << cells.size() << " != " << (int)flux_map->ElementSize(f) << "\n";
        Exceptions::amanzi_throw(msg);
      }

      double linear(0.0);
      std::vector<double> val(num_vec + 1, 0.0); // extra vector is for water exchange
      int pos = Operators::UniqueIndexFaceToCells(*mesh_out_, f, cells[0]);

      for (int j = 0; j != cells.size(); ++j) {
        auto [faces,dirs] = mesh_out_->getCellFacesAndDirections(cells[j]);

        for (int i = 0; i < faces.size(); i++) {
          if (f == faces[i]) {
            int g = flux_map->FirstPointInElement(f);
            double fln = flux[0][g + (pos + j) % 2] * dirs[i];

            if (fln >= 0) {
              linear += fln;
              for (int k = 0; k < num_vec; ++k) {
                // solute flux is water flux * concentration
                val[k] += field_out[k][cells[j]] * fln / mesh_->getCellVolume(c);
              }
            } else {
              val[num_vec] -= fln * (t1 - t0);
            }
            break;
          }
        }
      }
      value_[c] = val;
      linear_term_[c] = linear;
    }

  } else if (submodel_ == "flux exchange") {
    const auto& flux =
      *S_->Get<CompositeVector>(flux_key_, copy_flux_tag_).ViewComponent("face", true);
    const auto& flux_map =
      S_->Get<CompositeVector>(flux_key_, copy_flux_tag_).Map().Map("face", true);

    AmanziMesh::Entity_ID_List cells;

    // loop over cells on the manifold
    for (auto c : *entity_ids_) {
      AmanziMesh::Entity_ID f = mesh_->getEntityParent(AmanziMesh::Entity_kind::CELL, c);

      cells = mesh_out_->getFaceCells(f, AmanziMesh::Parallel_type::ALL);

      if (cells.size() != flux_map->ElementSize(f)) {
        msg << "Number of flux DOFs doesn't equal to the number of cells sharing a face: "
            << cells.size() << " != " << (int)flux_map->ElementSize(f) << "\n";
        Exceptions::amanzi_throw(msg);
      }

      std::vector<double> val(1, 0.0);
      int pos = Operators::UniqueIndexFaceToCells(*mesh_out_, f, cells[0]);

      for (int j = 0; j != cells.size(); ++j) {
        auto [faces, dirs] = mesh_out_->getCellFacesAndDirections(cells[j]);

        for (int i = 0; i < faces.size(); i++) {
          if (f == faces[i]) {
            int g = flux_map->FirstPointInElement(f);
            double fln = flux[0][g + (pos + j) % 2] * dirs[i];
            val[0] += fln / mesh_->getCellVolume(c);
            break;
          }
        }
      }
      value_[c] = val;
    }

  } else if (submodel_ == "field") {
    const auto& field_out =
      *S_->Get<CompositeVector>(field_out_key_, copy_field_out_tag_).ViewComponent("cell", true);

    int num_vec = field_out.NumVectors();
    std::vector<double> val(num_vec);

    // Loop over faces (owned + ghosted) in the space restricted to the manifold.
    // The set of these faces could be bigger then the set of manifold cells (owned + ghosted)
    for (auto f = entity_ids_->begin(); f != entity_ids_->end(); ++f) {
      auto it = reverse_map_.find(*f);
      if (it == reverse_map_.end()) continue;

      int c = it->second;
      for (int k = 0; k < num_vec; ++k) val[k] = field_out[k][c];
      value_[*f] = val;
    }

  } else if (submodel_ == "conserved quantity") {
    const auto& field_out =
      *S_->Get<CompositeVector>(field_out_key_, copy_field_out_tag_).ViewComponent("cell", true);
    const auto& field_cons =
      *S_->Get<CompositeVector>(field_cons_key_, copy_field_cons_tag_).ViewComponent("cell", true);

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
      AmanziMesh::Entity_ID_List cells;
      cells = mesh_->getFaceCells(f, AmanziMesh::Parallel_type::OWNED);
      AMANZI_ASSERT(cells.size() == 1);

      auto [faces, dirs] = mesh_->getCellFacesAndDirections(cells[0]);
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

} // namespace Amanzi

#endif
