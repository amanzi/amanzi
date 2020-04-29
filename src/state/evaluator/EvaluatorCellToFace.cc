/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/
//! Calculates a face value from cell values.

#include "errors.hh"
#include "EvaluatorCellToFace.hh"

namespace Amanzi {

EvaluatorCellToFace::EvaluatorCellToFace(Teuchos::ParameterList& plist) :
    EvaluatorSecondary(plist) {
  algorithm_ = plist.get<std::string>("averaging algorithm");
  if (algorithm_ != "harmonic" &&
      algorithm_ != "arithmetic" &&
      algorithm_ != "geometric") {
    Errors::Message msg;
    msg << "EvaluatorCellToFace: averaging algorithm \"" << algorithm_ << "\" not valid.  Valid are: \"harmonic\", \"arithmetic\", \"geometric\"";
    throw(msg);
  }
}

Evaluator&
EvaluatorCellToFace::operator=(const Evaluator& other)
{
  if (this != &other) {
    const EvaluatorCellToFace* other_p =
      dynamic_cast<const EvaluatorCellToFace*>(&other);
    AMANZI_ASSERT(other_p != NULL);
    *this = *other_p;
  }
  return *this;
}

void
EvaluatorCellToFace::EnsureCompatibility(State& S)
{
  AMANZI_ASSERT(my_keys_.size() == 1);
  AMANZI_ASSERT(dependencies_.size() == 1);
  auto my_key = my_keys_[0];
  auto& my_fac = S.Require<CompositeVector,CompositeVectorSpace>(my_key.first, my_key.second, my_key.first);

  if (my_fac.Mesh().get()) {
    my_fac.SetComponent("face", AmanziMesh::FACE, 1)->SetGhosted();
    auto& dep_fac = S.Require<CompositeVector,CompositeVectorSpace>(dependencies_[0].first, dependencies_[0].second);
    dep_fac.SetMesh(my_fac.Mesh());
    dep_fac.AddComponent("cell", AmanziMesh::CELL, 1)->SetGhosted();
    S.RequireEvaluator(dependencies_[0].first, dependencies_[0].second).EnsureCompatibility(S);
  }

  EnsureCompatibility_Flags_(S);
}

void
EvaluatorCellToFace::Update_(State& S)
{
  // hacked together implementation...
  const auto& cells = S.Get<CompositeVector>(dependencies_[0].first, dependencies_[0].second);
  cells.ScatterMasterToGhosted("cell");

  { // scope for views
    const AmanziMesh::Mesh* m = cells.getMap()->Mesh().get();
    auto cells_v = cells.ViewComponent("cell");
    auto faces_v = S.GetW<CompositeVector>(my_keys_[0].first, my_keys_[0].second,
            my_keys_[0].first).ViewComponent("face", false);

    if (algorithm_ == "harmonic") {
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
    } else if (algorithm_ == "arithmetic") {
      Kokkos::parallel_for(
          "EvaluatorCellToFace: arithmetic",
          faces_v.extent(0),
          KOKKOS_LAMBDA(const int& f) {
            AmanziMesh::Entity_ID_View cells;
            m->face_get_cells(f, AmanziMesh::Parallel_type::ALL, cells);
            faces_v(f,0) = cells.extent(0) == 1
                           ? cells_v(cells(0),0)
                           : (cells_v(cells(0),0) + cells_v(cells(1),0)) / 2.0;
          });

    } else if (algorithm_ == "geometric") {
      Kokkos::parallel_for(
          "EvaluatorCellToFace: geometric",
          faces_v.extent(0),
          KOKKOS_LAMBDA(const int& f) {
            AmanziMesh::Entity_ID_View cells;
            m->face_get_cells(f, AmanziMesh::Parallel_type::ALL, cells);
            faces_v(f,0) = cells.extent(0) == 1
                           ? cells_v(cells(0),0)
                           : sqrt(cells_v(cells(0),0) * cells_v(cells(1),0));
          });
    } else {
      AMANZI_ASSERT(false);
    }
  }
}

void
EvaluatorCellToFace::UpdateDerivative_(State& S, const Key& wrt_key, const Key& wrt_tag)
{
  AMANZI_ASSERT(false);
}


} // namespace Amanzi
