/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

//! Evaluates a dispersion tensor.
/*!

The dispersion tensor may be parameterized by a variety of different models,
and is specified by region to allow different models across, e.g. stream vs
sheet flow, or fracture vs matrix.

Note that this computes the dispersion coefficient tensor, units [m^2 s^-1].

It does NOT include the water content term which multiplies this tensor to
become the diffusion tensor.  This is different from Amanzi, which uses the MDM
class to include the water content term.

`"evaluator type`" = `"dispersion tensor`"

.. _evaluator-mdm-spec:
.. admonition: evaluator-mdm-spec

   * `"mechanical dispersion parameters`"
     ``[mechanical-dispersion-typed-spec-list]`` List of region-based
     models defining the dispersion model.

   * `"is surface`" ``[bool]`` **optional** Default is set if "surface" is in
     the domain name of this variable's key, e.g. in the case
     "surface-dispersion_coefficient".

   KEYS:

   - `"velocity`"
   - `"porosity`" Note this is only used if `"is surface`" = False.
   
*/

#pragma once

#include "dbc.hh"
#include "EvaluatorSecondary.hh"
#include "MDMPartition.hh"

namespace Amanzi {
namespace Transport {

class EvaluatorMDM : public EvaluatorSecondary
{
 public:
  // constructor format for all derived classes
  explicit EvaluatorMDM(Teuchos::ParameterList& plist);

  EvaluatorMDM(const EvaluatorMDM& other) = default;
  virtual Teuchos::RCP<Evaluator> Clone() const override {
    return Teuchos::rcp(new EvaluatorMDM(*this));
  }

 protected:
  // Required methods from EvaluatorSecondaryMonotypeCV
  virtual void Update_(State& S) override;

  virtual void UpdateDerivative_(State& S,
          const Key& wrt_key,
          const Tag& wrt_tag) override
  {
    AMANZI_ASSERT(false); // not implemented
  }

  virtual void EnsureCompatibility(State& S) override;

  virtual bool IsDifferentiableWRT(const State& S, const Key& wrt_key, const Tag& wrt_tag) const override
  {
    // derivatives not implemented
    return false;
  }

 protected:
  Teuchos::RCP<MDMPartition> mdms_;

  bool is_surface_;
  Key poro_key_;
  Key velocity_key_;

 private:
  static Utils::RegisteredFactory<Evaluator, EvaluatorMDM> reg_;


};

} // namespace Transport
} // namespace Amanzi
