/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov
      
*/
//! Calculates a face value from cell values using upwind. 

#include "errors.hh"
#include "EvaluatorCellToFaceUpwind.hh"

namespace Amanzi {

EvaluatorCellToFaceUpwind::EvaluatorCellToFaceUpwind(Teuchos::ParameterList& plist) :
    EvaluatorSecondary(plist) {
  flux_key_ = plist.get<std::string>("flux direction key");
}

Evaluator&
EvaluatorCellToFaceUpwind::operator=(const Evaluator& other)
{
  if (this != &other) {
    const auto* other_p = dynamic_cast<const EvaluatorCellToFaceUpwind*>(&other);
    AMANZI_ASSERT(other_p != NULL);
    *this = *other_p;
  }
  return *this;
}

void
EvaluatorCellToFaceUpwind::EnsureCompatibility(State& S)
{
  AMANZI_ASSERT(my_keys_.size() == 1);
  AMANZI_ASSERT(dependencies_.size() == 1);

  auto my_key = my_keys_[0];
  auto& my_fac = S.Require<CompositeVector,CompositeVectorSpace>(my_key.first, my_key.second, my_key.first);

  if (my_fac.Mesh().get()) {
    my_fac.SetComponent("face", AmanziMesh::FACE, 1)->SetGhosted();

    for (auto deps : dependencies_) {
      auto& dep_fac = S.Require<CompositeVector,CompositeVectorSpace>(deps.first, deps.second);
      dep_fac.SetMesh(my_fac.Mesh());
      dep_fac.AddComponent("cell", AmanziMesh::CELL, 1)->SetGhosted();

      S.RequireEvaluator(deps.first, deps.second).EnsureCompatibility(S);

      bool has_derivs = false;
      for (auto keytag : my_keys_)
        has_derivs |= S.HasDerivativeSet(keytag.first, keytag.second);
      if (has_derivs) {
        for (const auto& deriv : S.GetDerivativeSet(my_key.first, my_key.second)) {
          auto wrt = Keys::splitKeyTag(deriv.first);
          S.RequireDerivative<CompositeVector, CompositeVectorSpace>(
              my_key.first, my_key.second, wrt.first, wrt.second, my_key.first).Update(my_fac);

          S.RequireDerivative<CompositeVector, CompositeVectorSpace>(
              deps.first, deps.second, wrt.first, wrt.second).Update(dep_fac);
        }
      }    
    }
  }

  EnsureCompatibility_Flags_(S);
}

void
EvaluatorCellToFaceUpwind::Update_(State& S)
{
  // hacked together implementation...
  const auto& cells = S.Get<CompositeVector>(dependencies_[0].first, dependencies_[0].second);
  cells.ScatterMasterToGhosted("cell");

  { // scope for views
    const AmanziMesh::Mesh* m = cells.getMap()->Mesh().get();
    auto cells_v = cells.ViewComponent("cell");
    auto faces_v = S.GetW<CompositeVector>(my_keys_[0].first, my_keys_[0].second,
            my_keys_[0].first).ViewComponent("face", false);

    const auto& flux = S.Get<CompositeVector>(flux_key_, "next").ViewComponent("face");

    Kokkos::parallel_for(
        "EvaluatorCellToFace: upwind",
        faces_v.extent(0),
        KOKKOS_LAMBDA(const int& f) {
          int dir;
          AmanziMesh::Entity_ID_View cells;
          m->face_get_cells(f, AmanziMesh::Parallel_type::ALL, cells);
          const auto& normal = m->face_normal(f, false, cells(0), &dir);

          faces_v(f,0) = (cells.extent(0) == 1)
                         ? cells_v(cells(0),0)
                         : (flux(f,0) * dir > 0) ? cells_v(cells(0),0) : cells_v(cells(1),0);
        });
  }
}

void
EvaluatorCellToFaceUpwind::UpdateDerivative_(State& S, const Key& wrt_key, const Key& wrt_tag)
{
  // hacked together implementation...
  const auto& cells = S.GetDerivative<CompositeVector>(dependencies_[0].first, dependencies_[0].second,
          wrt_key, wrt_tag);
  cells.ScatterMasterToGhosted("cell");

  { // scope for views
    const AmanziMesh::Mesh* m = cells.getMap()->Mesh().get();
    auto cells_v = cells.ViewComponent("cell");
    auto faces_v = S.GetDerivativeW<CompositeVector>(my_keys_[0].first, my_keys_[0].second,
            wrt_key, wrt_tag, my_keys_[0].first).ViewComponent("face", false);

    Kokkos::parallel_for(
        "EvaluatorCellToFace: harmonic",
        faces_v.extent(0),
        KOKKOS_LAMBDA(const int& f) {
          AmanziMesh::Entity_ID_View cells;
          m->face_get_cells(f, AmanziMesh::Parallel_type::ALL, cells);
          faces_v(f,0) = cells.extent(0) == 1
                         ? cells_v(cells(0),0)
                         : 1.0 / (1./cells_v(cells(0),0) + 1./cells_v(cells(1),0));
        });
  }
}

} // namespace Amanzi
