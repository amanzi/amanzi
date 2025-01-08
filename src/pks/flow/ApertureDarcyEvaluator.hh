/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
           Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Flow PK

  The linearized fracture aperture evaluator for Darcy model.
*/

#ifndef AMANZI_APERTURE_EVALUATOR_DARCY_HH_
#define AMANZI_APERTURE_EVALUATOR_DARCY_HH_

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"
#include "Key.hh"
#include "Tag.hh"

#include "ApertureModel.hh"
#include "ApertureModelPartition.hh"

namespace Amanzi {
namespace Flow {

class ApertureDarcyEvaluator
  : public EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace> {
 public:
  // constructor format for all derived classes
  explicit ApertureDarcyEvaluator(Teuchos::ParameterList& plist);
  ApertureDarcyEvaluator(const ApertureDarcyEvaluator& other);

  // required inteface functions
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& results) override;

  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& results) override;

 protected:
  virtual void EnsureCompatibility_ToDeps_(State& S) override{};

 protected:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  Key ref_pressure_key_, ref_aperture_key_, pressure_key_, compliance_key_;

 private:
  static Utils::RegisteredFactory<Evaluator, ApertureDarcyEvaluator> reg_;
};

} // namespace Flow
} // namespace Amanzi

#endif
