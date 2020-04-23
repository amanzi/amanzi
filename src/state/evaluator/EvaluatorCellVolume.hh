/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:

*/

//!

/*
  State

  License: BSD
  Author: Ethan Coon

*/

#ifndef AMANZI_INDEPENDENT_CELL_VOLUME_EVALUATOR_HH_
#define AMANZI_INDEPENDENT_CELL_VOLUME_EVALUATOR_HH_

#include "EvaluatorIndependent.hh"

namespace Amanzi {

class EvaluatorCellVolume
  : public EvaluatorIndependent<CompositeVector, CompositeVectorSpace> {
 public:
  // ---------------------------------------------------------------------------
  // Constructors, assignement operators, etc
  // ---------------------------------------------------------------------------
  using EvaluatorIndependent<CompositeVector,
                             CompositeVectorSpace>::EvaluatorIndependent;
  EvaluatorCellVolume(const EvaluatorCellVolume& other) = default;

  Teuchos::RCP<Evaluator> Clone() const override
  {
    return Teuchos::rcp(new EvaluatorCellVolume(*this));
  }

  virtual std::string name() const override { return "cell volume"; }
  
 protected:
  // ---------------------------------------------------------------------------
  // Does the actual work to update the value in the state.
  // ---------------------------------------------------------------------------
  virtual void Update_(State& S) override;

 private:
  static Utils::RegisteredFactory<Evaluator, EvaluatorCellVolume> fac_;
};

//
// Provide limited scope to private data
//
namespace Impl {
  template<class View_type>
  void copyCellVolume(const AmanziMesh::Mesh* mesh, View_type& v) {
    Kokkos::parallel_for("EvaluatorCellVolume copy",
                         v.extent(0),
                         KOKKOS_LAMBDA(const int& c) {
                           v(c,0) = mesh->cell_volume(c);
                         });
  }

  template<class View_type>
  void copyFaceArea(const AmanziMesh::Mesh* mesh, View_type& v) {
    Kokkos::parallel_for("EvaluatorCellVolume copy",
                         v.extent(0),
                         KOKKOS_LAMBDA(const int& c) {
                           v(c,0) = mesh->face_area(c);
                         });
  }
} // namespace Impl
  
  
} // namespace Amanzi

#endif
