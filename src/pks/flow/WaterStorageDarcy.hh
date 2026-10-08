/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Flow PK

  Field evaluator for water storage which is the conserved quantity
  in the Darcy equation.

  Wrapping this conserved quantity water storage (WS) as a field evaluator
  makes it easier to take derivatives, keep updated, and the like.
  The equation for this is simply:

    WS = a * Ss * p / g

  where a is the optional fracture aperture (default is a = 1), Ss is the
  specific storage, p is the liquid pressure, and g is the gravity magnitude.
*/


#ifndef AMANZI_FLOW_WATER_STORAGE_DARCY_EVALUATOR_HH_
#define AMANZI_FLOW_WATER_STORAGE_DARCY_EVALUATOR_HH_

#include "Teuchos_ParameterList.hpp"

#include "EvaluatorSecondaryMonotype.hh"
#include "Factory.hh"
#include "Tag.hh"

namespace Amanzi {
namespace Flow {

class WaterStorageDarcy : public EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace> {
 public:
  explicit WaterStorageDarcy(Teuchos::ParameterList& plist);
  WaterStorageDarcy(const WaterStorageDarcy& other);

  // required inteface functions
  virtual Teuchos::RCP<Evaluator> Clone() const override;

  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& results) override;

  virtual void EvaluatePartialDerivative_(const State& S,
                                          const Key& wrt_key,
                                          const Tag& wrt_tag,
                                          const std::vector<CompositeVector*>& results) override;

 protected:
  bool aperture_;
  Key specific_storage_key_, aperture_key_, pressure_key_;

 private:
  static Utils::RegisteredFactory<Evaluator, WaterStorageDarcy> reg_;
};

} // namespace Flow
} // namespace Amanzi

#endif
