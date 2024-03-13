/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@ornl.gov)
*/

//! Adapts EvaluatorModelCVByMaterial to evaluate multi-dof things.

#pragma once

#include "Factory.hh"
#include "EvaluatorModelCVByMaterial.hh"

namespace Amanzi {

template <template <class, class> class Model, class Device_type = DefaultDevice>
class EvaluatorMultiDOFModelCVByMaterial
  : public EvaluatorModelCVByMaterial<Model, Device_type> {
 public:
  // constructor via plist
  using EvaluatorModelCVByMaterial<Model, Device_type>::EvaluatorModelCVByMaterial;

  // copy constructor
  EvaluatorMultiDOFModelCVByMaterial(const EvaluatorMultiDOFModelCVByMaterial& other) = default;

  // virtual copy constructor
  virtual Teuchos::RCP<Evaluator> Clone() const override {
    return Teuchos::rcp(new EvaluatorMultiDOFModelCVByMaterial(*this));
  }

 protected:
  // custom EC used to set subfield names and number of dofs
  virtual void EnsureCompatibility_Structure_(State& S) override {
    for (const auto& key_tag : my_keys_) {
      S.Require<CompositeVector, CompositeVectorSpace>(key_tag.first, key_tag.second)
        .SetMesh(S.GetMesh(Keys::getDomain(key_tag.first)))
        ->SetGhosted()
        ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, Model_type::n_dofs);
      S.GetRecordSetW(key_tag.first).set_subfieldnames(Model_type::subfield_names);
    }
  }

  // custom EC used because deps have 1 component not N
  virtual void EnsureCompatibility_ToDeps_(State& S) override {
    for (const auto& dep : dependencies_) {
      S.Require<CompositeVector, CompositeVectorSpace>(dep.first, dep.second)
        .SetMesh(S.GetMesh(Keys::getDomain(dep.first)))
        ->SetGhosted()
        ->SetComponent("cell", AmanziMesh::Entity_kind::CELL, 1);
    }
  }

 private:
  using EvaluatorModelCVByMaterial<Model, Device_type>::eval_type;
  using EvaluatorModelCVByMaterial<Model, Device_type>::my_keys_;
  using EvaluatorModelCVByMaterial<Model, Device_type>::dependencies_;
  using typename EvaluatorModelCVByMaterial<Model, Device_type>::Model_type;
  static Utils::RegisteredFactory<Evaluator, EvaluatorMultiDOFModelCVByMaterial> reg_;
};

} // namespace Amanzi
