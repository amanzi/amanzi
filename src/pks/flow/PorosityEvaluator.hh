/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
           Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Flow PK

  The porosity evaluator simply calls the porosity model with
  the correct arguments.
*/

#ifndef AMANZI_FLOW_POROSITY_EVALUATOR_HH_
#define AMANZI_FLOW_POROSITY_EVALUATOR_HH_

#include <utility>

#include "Factory.hh"
#include "EvaluatorSecondaryMonotype.hh"
#include "Tag.hh"

#include "Porosity.hh"
#include "PorosityModelPartition.hh"

namespace Amanzi {
namespace Flow {

class PorosityEvaluator : public EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace> {
 public:
  // constructor format for all derived classes
  explicit PorosityEvaluator(Teuchos::ParameterList& plist,
                             Teuchos::RCP<PorosityModelPartition> pom);
  PorosityEvaluator(const PorosityEvaluator& other);

  // required inteface functions
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& results) override;

  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& results) override;

  std::pair<double, double> getThermalCoefficients(int c)
  {
    return pom_->second[(*pom_->first)[c]]->getThermalCoefficients();
  }

 protected:
  void InitializeFromPlist_();

 protected:
  Teuchos::RCP<const AmanziMesh::Mesh> mesh_;
  Teuchos::RCP<PorosityModelPartition> pom_;
  Key pressure_key_, temperature_key_, strain_key_, biot_key_;

  bool poroelasticity_, thermoelasticity_;

 private:
  static Utils::RegisteredFactory<Evaluator, PorosityEvaluator> reg_;
};

} // namespace Flow
} // namespace Amanzi

#endif
