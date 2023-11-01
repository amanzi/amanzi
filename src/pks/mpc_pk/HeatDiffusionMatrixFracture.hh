/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*!

Heat conduction coefficient between fracture and matrix:

  D_fm = D_m / d_mf

where D_m is matrix heat thermal conductivity, and d_mf is distance
between fracture and matrix cell centroids. This is a FV approximation
of the heat flux.
*/


#ifndef AMANZI_HEAT_DIFFUSION_MATRIX_FRACTURE_HH_
#define AMANZI_HEAT_DIFFUSION_MATRIX_FRACTURE_HH_

#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"

#include "EvaluatorSecondary.hh"

namespace Amanzi {

class HeatDiffusionMatrixFracture : public EvaluatorSecondary {
 public:
  HeatDiffusionMatrixFracture(Teuchos::ParameterList& plist);
  HeatDiffusionMatrixFracture(const HeatDiffusionMatrixFracture& other) = default;

  virtual Teuchos::RCP<Evaluator> Clone() const override
  {
    return Teuchos::rcp(new HeatDiffusionMatrixFracture(*this));
  }

  virtual void EnsureCompatibility(State& S) override;

 protected:
  virtual void Update_(State& S) override;
  virtual void UpdateDerivative_(State& S, const Key& wrt_key, const Tag& wrt_tag) override
  {
    AMANZI_ASSERT(false);
  }

 private:
  Key domain_;
  Key conductivity_key_, aperture_key_;

  static Utils::RegisteredFactory<Evaluator, HeatDiffusionMatrixFracture> reg_;
};

} // namespace Amanzi

#endif
