/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*!

Molecular diffusion coefficeint between fracture and matrix.

  D_fm = s_m tau_m phi_m * D_m / d_mf

where D_m is matrix diffusion coefficient, s_m is saturation,
tau_m is tortuosity, and d_mf is distance between fracture and
matrix cell centroids. This is a FV approximation of the solute
flux.
*/


#ifndef AMANZI_SOLUTE_DIFFUSION_MATRIX_FRACTURE_HH_
#define AMANZI_SOLUTE_DIFFUSION_MATRIX_FRACTURE_HH_

#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"

#include "EvaluatorSecondary.hh"


namespace Amanzi {

class SoluteDiffusionMatrixFracture : public EvaluatorSecondary {
 public:
  SoluteDiffusionMatrixFracture(Teuchos::ParameterList& plist);
  SoluteDiffusionMatrixFracture(const SoluteDiffusionMatrixFracture& other) = default;

  virtual Teuchos::RCP<Evaluator> Clone() const override
  {
    return Teuchos::rcp(new SoluteDiffusionMatrixFracture(*this));
  }

  virtual void EnsureCompatibility(State& S) override;

  // modifiers
  // -- WIP for subvectors
  void set_subvector(double mol_diff)
  {
    mol_diff_ = mol_diff;
    requests_.clear();
  }

 protected:
  virtual void Update_(State& S) override;
  virtual void UpdateDerivative_(State& S, const Key& wrt_key, const Tag& wrt_tag) override
  {
    AMANZI_ASSERT(false);
  }

 private:
  Key domain_;
  Key saturation_key_, tortuosity_key_, porosity_key_, aperture_key_;
  double mol_diff_;

  static Utils::RegisteredFactory<Evaluator, SoluteDiffusionMatrixFracture> reg_;
};

} // namespace Amanzi

#endif
