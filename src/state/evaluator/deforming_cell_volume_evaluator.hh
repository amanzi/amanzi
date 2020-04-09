/*
  Copyright 2010-201x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
      Ethan Coon
*/

//!

#ifndef AMANZI_DEFORMING_CELL_VOLUME_FIELD_EVALUATOR_
#define AMANZI_DEFORMING_CELL_VOLUME_FIELD_EVALUATOR_

#include "EvaluatorSecondaryMonotype.hh"
#include "Evaluator_Factory.hh"

namespace Amanzi {

class DeformingCellVolumeEvaluator : public EvaluatorSecondaryMonotype {
 public:
  // constructor format for all derived classes
  explicit DeformingCellVolumeEvaluator(Teuchos::ParameterList& plist);

  // copy constructors
  DeformingCellVolumeEvaluator(const DeformingCellVolumeEvaluator& other);
  virtual Teuchos::RCP<Evaluator> Clone() const;

  // setup of data and class
  virtual void EnsureCompatibility(const Teuchos::Ptr<State>& S);

 protected:
  // Required methods from EvaluatorSecondaryMonotype
  virtual void EvaluateField_(const Teuchos::Ptr<State>& S,
                              const Teuchos::Ptr<CompositeVector>& result);
  virtual void
  EvaluateFieldPartialDerivative_(const Teuchos::Ptr<State>& S, Key wrt_key,
                                  const Teuchos::Ptr<CompositeVector>& result);

 protected:
  Key my_mesh_;

 private:
  static Utils::RegisteredFactory<Evaluator, DeformingCellVolumeEvaluator> fac_;
};

} // namespace Amanzi

#endif
