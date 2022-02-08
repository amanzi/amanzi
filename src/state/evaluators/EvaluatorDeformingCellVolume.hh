/*
  State

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon

  A field evaluator for an changing cell volume.
*/

#ifndef AMANZI_STATE_DEFORMING_CELL_VOLUME_FIELD_EVALUATOR_
#define AMANZI_STATE_DEFORMING_CELL_VOLUME_FIELD_EVALUATOR_

#include "Evaluator_Factory.hh"
#include "EvaluatorSecondary.hh"

namespace Amanzi {

class EvaluatorDeformingCellVolume : public EvaluatorSecondary {
 public:
  // constructor format for all derived classes
  explicit EvaluatorDeformingCellVolume(Teuchos::ParameterList& plist);

  // copy constructors
  EvaluatorDeformingCellVolume(const EvaluatorDeformingCellVolume &other) = default;
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  // not differentiable WRT to anything
  virtual bool IsDifferentiableWRT(const State& S, const Key& wrt_key,
          const Tag& wrt_tag) const override { return false; }

  virtual void EnsureCompatibility(State& S) override;

 protected:

  virtual void Update_(State &S) override;
  virtual void UpdateDerivative_(State& S, const Key& wrt_key, const Tag& wrt_tag) override {};

 protected:
  Key my_mesh_, my_key_;

 private:
  static Utils::RegisteredFactory<Evaluator, EvaluatorDeformingCellVolume> fac_;
};

} // namespace Amanzi

#endif
