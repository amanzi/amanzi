/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

/*!

An evaluator that depends upon a static mesh, and evaluates various functions
for quantities in or derived from the mesh.  Note these _could_ stay in the
mesh and not go into a vector in many cases (especially, for instance, cell
volume) but for a combination of laziness (easier to write
Vector::elementWiseMultiply() than manually write a loop), generality (there
are times when some quantities could come from a mesh or could come from
another way), or documentation (have cell volumes in vis is really useful),
these have been put into vectors and saved in state instead of being
recalculated from the mesh each time.

*/

#pragma once

#include "MeshDefs.hh"
#include "EvaluatorSecondary.hh"
#include "Evaluator_Factory.hh"

namespace Amanzi {

template<class Function_type>
class EvaluatorSecondaryMeshedQuantity
  : public EvaluatorSecondary {
 public:
  // ---------------------------------------------------------------------------
  // Constructors, assignement operators, etc
  // ---------------------------------------------------------------------------
  EvaluatorSecondaryMeshedQuantity(Teuchos::ParameterList& plist)
      : EvaluatorSecondary(plist)
  {}

  EvaluatorSecondaryMeshedQuantity(const EvaluatorSecondaryMeshedQuantity& other) = default;

  virtual Teuchos::RCP<Evaluator> Clone() const override
  {
    return Teuchos::rcp(new EvaluatorSecondaryMeshedQuantity(*this));
  }

  virtual std::string name() const override { return "meshed quantity evaluator"; }

  virtual void EnsureCompatibility(State& S) override {
    AMANZI_ASSERT(my_keys_.size() == 1);
    std::string domain = Keys::getDomain(my_keys_[0].first);

    S.Require<CompositeVector,CompositeVectorSpace>(my_keys_[0].first,
            my_keys_[0].second, my_keys_[0].first)
        .SetMesh(S.GetMesh(domain))
        ->SetComponent(AmanziMesh::entity_kind_string(Function_type::component),
                      Function_type::component, 1)
        ->SetGhosted(true);
    S.RequireEvaluator(Keys::getKey(domain, "mesh"));
  }
  
 protected:
  // ---------------------------------------------------------------------------
  // Does the actual work to update the value in the state.
  // ---------------------------------------------------------------------------
  virtual void Update_(State& S) override {
    std::string domain = Keys::getDomain(my_keys_[0].first);
    auto mesh = S.GetMesh(domain);
    auto v = S.GetPtrW<CompositeVector>(my_keys_[0].first, my_keys_[0].second, my_keys_[0].first);
    Function_type f(mesh.get(), v);
    f.Compute();
  }

  virtual void UpdateDerivative_(State& S, const Key& wrt_key, const Key& wrt_tag) override
  {
    AMANZI_ASSERT(false);
  }
  
 private:
  static Utils::RegisteredFactory<Evaluator, EvaluatorSecondaryMeshedQuantity<Function_type>> fac_;
};

//
// Provide limited scope to private data
//
namespace Impl {
  struct CellVolume {
    CellVolume(const AmanziMesh::Mesh* mesh_,
               Teuchos::RCP<CompositeVector>& v_)
        : mesh(mesh_),
          v(v_) {}

    void Compute() {
      auto vv = v->ViewComponent("cell", false);
      Kokkos::parallel_for(
          "EvaluatorCellVolume",
          vv.extent(0),
          KOKKOS_LAMBDA(const int& c) {
            vv(c,0) = mesh->cell_volume(c);
          });
    }

    const static AmanziMesh::Entity_kind component = AmanziMesh::Entity_kind::CELL;

    const AmanziMesh::Mesh* mesh;
    Teuchos::RCP<CompositeVector> v;
  };


  struct FaceArea {
    FaceArea(const AmanziMesh::Mesh* mesh_,
               Teuchos::RCP<CompositeVector>& v_)
        : mesh(mesh_),
          v(v_) {}

    void Compute() {
      auto vv = v->ViewComponent("face", false);
      Kokkos::parallel_for(
          "EvaluatorFaceArea",
          vv.extent(0),
          KOKKOS_LAMBDA(const int& c) {
            vv(c,0) = mesh->face_area(c);
          });
    }

    const static AmanziMesh::Entity_kind component = AmanziMesh::Entity_kind::FACE;
    const AmanziMesh::Mesh* mesh;
    Teuchos::RCP<CompositeVector> v;
  };

  struct CellElevation {
    CellElevation(const AmanziMesh::Mesh* mesh_,
               Teuchos::RCP<CompositeVector>& v_)
        : mesh(mesh_),
          v(v_) {}

    void Compute() {
      auto vv = v->ViewComponent("cell", false);
      const AmanziMesh::Mesh* parent = mesh->parent().get();
      int d = parent->space_dimension()-1;

      if (parent->space_dimension() == mesh->space_dimension()) {
        Kokkos::parallel_for(
            "EvaluatorCellElevation",
            vv.extent(0),
            KOKKOS_LAMBDA(const int& c) {
              vv(c,0) = parent->cell_centroid(mesh->entity_get_parent(AmanziMesh::CELL,c))[d];
            });
      } else if (parent->space_dimension() == mesh->space_dimension() + 1) {
        Kokkos::parallel_for(
            "EvaluatorCellElevation",
            vv.extent(0),
            KOKKOS_LAMBDA(const int& c) {
              vv(c,0) = parent->face_centroid(mesh->entity_get_parent(AmanziMesh::CELL,c))[d];
            });
      } else {
        AMANZI_ASSERT(false);
      }        
    }

    const static AmanziMesh::Entity_kind component = AmanziMesh::Entity_kind::CELL;
    const AmanziMesh::Mesh* mesh;
    Teuchos::RCP<CompositeVector> v;
  };

  struct SlopeMagnitude {
    SlopeMagnitude(const AmanziMesh::Mesh* mesh_,
               Teuchos::RCP<CompositeVector>& v_)
        : mesh(mesh_),
          v(v_) {}

    void Compute() {
      auto vv = v->ViewComponent("cell", false);
      const AmanziMesh::Mesh* parent = mesh->parent().get();

      if (parent->space_dimension() == mesh->space_dimension()) {
        Kokkos::parallel_for(
            "EvaluatorCellElevation",
            vv.extent(0),
            KOKKOS_LAMBDA(const int& c) {
              // Now slope.
              auto f = mesh->entity_get_parent(AmanziMesh::CELL,c);
              AmanziGeometry::Point n = parent->face_normal(f);

              // -- S = || n - (n dot z) z || / | n dot z |
              vv(c,0) = sqrt(pow(n[0],2) + pow(n[1],2)) / fabs(n[2]);

            });
      } else if (parent->space_dimension() == mesh->space_dimension() + 1) {
        throw("not yet implemented");
      } else {
        AMANZI_ASSERT(false);
      }        
    }

    const static AmanziMesh::Entity_kind component = AmanziMesh::Entity_kind::CELL;
    const AmanziMesh::Mesh* mesh;
    Teuchos::RCP<CompositeVector> v;
  };


} // namespace Impl

using EvaluatorCellVolume = EvaluatorSecondaryMeshedQuantity<Impl::CellVolume>;
using EvaluatorFaceArea = EvaluatorSecondaryMeshedQuantity<Impl::FaceArea>;
using EvaluatorMeshElevation = EvaluatorSecondaryMeshedQuantity<Impl::CellElevation>;
using EvaluatorMeshSlopeMagnitude = EvaluatorSecondaryMeshedQuantity<Impl::SlopeMagnitude>;



} // namespace Amanzi


