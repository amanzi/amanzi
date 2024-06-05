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

template <template <AmanziMesh::Entity_kind> class Function_type>
class EvaluatorSecondaryMeshedQuantity : public EvaluatorSecondary {
 public:
  // ---------------------------------------------------------------------------
  // Constructors, assignement operators, etc
  // ---------------------------------------------------------------------------
  EvaluatorSecondaryMeshedQuantity(const Teuchos::RCP<Teuchos::ParameterList>& plist)
    : EvaluatorSecondary(plist), inited_(false)
  {}

  EvaluatorSecondaryMeshedQuantity(const EvaluatorSecondaryMeshedQuantity& other) = default;

  virtual Teuchos::RCP<Evaluator> Clone() const override
  {
    return Teuchos::rcp(new EvaluatorSecondaryMeshedQuantity(*this));
  }

  static const std::string eval_type;
  virtual std::string getType() const override { return eval_type; }

  virtual void EnsureCompatibility(State& S) override
  {
    if (!inited_) {
      AMANZI_ASSERT(my_keys_.size() == 1);
      domain_ = Keys::getDomain(my_keys_[0].first);
      domain3d_ = domain_;
      if (S.HasMesh(domain_ + "_3d")) domain3d_ = domain_ + "_3d";
      if (S.HasMesh(domain_ + "_3D")) domain3d_ = domain_ + "_3D";

      S.Require<CompositeVector, CompositeVectorSpace>(
         my_keys_[0].first, my_keys_[0].second, my_keys_[0].first)
        .SetMesh(S.GetMesh(domain_));

      S.RequireEvaluator(Keys::getKey(domain_, "mesh"), my_keys_[0].second);
      inited_ = true;
    }

    EnsureCompatibility_Flags_(S);
  }

 protected:
  // ---------------------------------------------------------------------------
  // Does the actual work to update the value in the state.
  // ---------------------------------------------------------------------------
  virtual void Update_(State& S) override
  {
    auto mesh = Function_type<AmanziMesh::Entity_kind::CELL>::needs_3d ? S.GetMesh(domain3d_) :
                                                                         S.GetMesh(domain_);
    auto v = S.GetPtrW<CompositeVector>(my_keys_[0].first, my_keys_[0].second, my_keys_[0].first);
    for (const auto& comp : *v) {
      if (comp == "cell") {
        Function_type<AmanziMesh::Entity_kind::CELL> f(*mesh, v);
        f.Compute();
      } else if (comp == "face") {
        Function_type<AmanziMesh::Entity_kind::FACE> f(*mesh, v);
        f.Compute();
      } else if (comp == "boundary_face") {
        Function_type<AmanziMesh::Entity_kind::BOUNDARY_FACE> f(*mesh, v);
        f.Compute();
      } else if (comp == "edge") {
        Function_type<AmanziMesh::Entity_kind::EDGE> f(*mesh, v);
        f.Compute();
      }
    }
    v->scatterMasterToGhosted();
  }

  virtual void UpdateDerivative_(State& S, const Key& wrt_key, const Tag& wrt_tag) override
  {
    AMANZI_ASSERT(false);
  }

 protected:
  bool inited_;
  Key domain_;
  Key domain3d_;

 private:
  static Utils::RegisteredFactory<Evaluator, EvaluatorSecondaryMeshedQuantity<Function_type>> reg_;
};

//
// Provide limited scope to private data
//
namespace Impl {

template <AmanziMesh::Entity_kind EK>
struct Extent {
  static const AmanziMesh::Entity_kind EntityKind = EK;
  static const bool needs_3d = false;

  Extent(const AmanziMesh::Mesh& mesh_, Teuchos::RCP<CompositeVector>& v_) : mesh(mesh_), v(v_) {}

  void Compute()
  {
    std::string compname = AmanziMesh::to_string(EntityKind);
    auto vv = v->viewComponent(compname, false);
    const auto& m = mesh;
    Kokkos::parallel_for(
      "EvaluatorMeshedEntityExtent", vv.extent(0), KOKKOS_LAMBDA(const int c) {
        vv(c, 0) = m.getExtent<EntityKind>(c);
      });
  }

  const AmanziMesh::Mesh& mesh;
  Teuchos::RCP<CompositeVector> v;
};


template <AmanziMesh::Entity_kind EK>
struct Elevation {
  static const AmanziMesh::Entity_kind EntityKind = EK;
  static const bool needs_3d = true;

  Elevation(const AmanziMesh::Mesh& mesh_, Teuchos::RCP<CompositeVector>& v_) : mesh(mesh_), v(v_)
  {}

  void Compute()
  {
    std::string compname = AmanziMesh::to_string(EntityKind);
    auto vv = v->viewComponent(compname, false);
    const auto& m = mesh;
    int d = mesh.getSpaceDimension() - 1;

    Kokkos::parallel_for(
      "EvaluatorMeshedEntityElevation", vv.extent(0), KOKKOS_LAMBDA(const int c) {
        vv(c, 0) = m.getCentroid<EntityKind>(c)[d];
      });
  }

  const AmanziMesh::Mesh& mesh;
  Teuchos::RCP<CompositeVector> v;
};


template <AmanziMesh::Entity_kind EK>
struct SlopeMagnitude {
  // only CELL is valid, but we need generality so we template it anyway
  static const AmanziMesh::Entity_kind EntityKind = EK;
  static const bool needs_3d = false;

  SlopeMagnitude(const AmanziMesh::Mesh& mesh_, Teuchos::RCP<CompositeVector>& v_)
    : mesh(mesh_), v(v_)
  {}

  void Compute()
  {
    AMANZI_ASSERT(EntityKind == AmanziMesh::Entity_kind::CELL);
    AMANZI_ASSERT(mesh.getParentMesh() != Teuchos::null);

    const AmanziMesh::Mesh& parent = *mesh.getParentMesh();
    auto vv = v->viewComponent("cell", false);
    const auto& m = mesh;

    Kokkos::parallel_for(
      "EvaluatorMeshedCellSlopeMagnitude", vv.extent(0), KOKKOS_LAMBDA(const int& c) {
        // face normal
        auto f = m.getEntityParent(AmanziMesh::Entity_kind::CELL, c);
        AmanziGeometry::Point n = parent.getFaceNormal(f);

        // -- S = || n - (n dot z) z || / | n dot z |
        vv(c, 0) = sqrt(pow(n[0], 2) + pow(n[1], 2)) / fabs(n[2]);
      });
  }

  const AmanziMesh::Mesh& mesh;
  Teuchos::RCP<CompositeVector> v;
};


template <AmanziMesh::Entity_kind EK>
struct Aspect {
  static const AmanziMesh::Entity_kind EntityKind = EK;
  static const bool needs_3d = false;

  Aspect(const AmanziMesh::Mesh& mesh_, Teuchos::RCP<CompositeVector>& v_) : mesh(mesh_), v(v_) {}

  void Compute()
  {
    AMANZI_ASSERT(EntityKind == AmanziMesh::Entity_kind::CELL);
    AMANZI_ASSERT(mesh.getParentMesh() != Teuchos::null);

    const AmanziMesh::Mesh& parent = *mesh.getParentMesh();
    auto vv = v->viewComponent("cell", false);
    const auto& m = mesh;

    Kokkos::parallel_for(
      "EvaluatorMeshedCellAspect", vv.extent(0), KOKKOS_LAMBDA(const int& c) {
        // compute the normal
        auto f = m.getEntityParent(AmanziMesh::Entity_kind::CELL, c);
        AmanziGeometry::Point n = parent.getFaceNormal(f);

        // and aspect
        if (n[0] > 0.0) {
          // right half
          if (n[1] > 0.0) {
            // upper right quadrant
            vv(c, 0) = atan(n[0] / n[1]);
          } else if (n[1] < 0.0) {
            // lower right quadrant
            vv(c, 0) = M_PI - atan(n[0] / -n[1]);
          } else {
            // due east
            vv(c, 0) = M_PI_2;
          }
        } else if (n[0] < 0.0) {
          // left half
          if (n[1] > 0.0) {
            // upper left quadrant
            vv(c, 0) = 2 * M_PI - atan(-n[0] / n[1]);
          } else if (n[1] < 0.0) {
            // lower left quadrant
            vv(c, 0) = M_PI + atan(n[0] / -n[1]);
          } else {
            // due west
            vv(c, 0) = 3 * M_PI_2;
          }
        } else {
          // north or south
          if (n[1] > 0.0) {
            vv(c, 0) = 0.0;
          } else if (n[1] < 0.0) {
            vv(c, 0) = M_PI;
          } else {
            // n is (0,0,1)
            vv(c, 0) = 0.;
          }
        }
      });
  }

  const AmanziMesh::Mesh& mesh;
  Teuchos::RCP<CompositeVector> v;
};


} // namespace Impl

using EvaluatorCellVolume = EvaluatorSecondaryMeshedQuantity<Impl::Extent>;
template <>
const std::string EvaluatorCellVolume::eval_type = "cell volume";

using EvaluatorMeshElevation = EvaluatorSecondaryMeshedQuantity<Impl::Elevation>;
template <>
const std::string EvaluatorMeshElevation::eval_type = "elevation";

using EvaluatorMeshSlopeMagnitude = EvaluatorSecondaryMeshedQuantity<Impl::SlopeMagnitude>;
template <>
const std::string EvaluatorMeshSlopeMagnitude::eval_type = "slope magnitude";

using EvaluatorMeshAspect = EvaluatorSecondaryMeshedQuantity<Impl::Aspect>;
template <>
const std::string EvaluatorMeshAspect::eval_type = "aspect";


} // namespace Amanzi
