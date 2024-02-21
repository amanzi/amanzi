/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

/*!

An "evaluator" that simply presents a vector as a patch.  This is useful for
things like boundary condition aggregators, which expect BCs to come in as
patches, must be given MultiVectors.  Two examples of this include (1) when
using a surface-subsurface flux vector as a Neuamnn BC for the subsurface in a
coupled surface-subsurface problem, and (2) when moving a field on one mesh to
provide a source/sink on another mesh, e.g. the subgrid hyporheic exchange flux
source to the surface system.

Note this does not really store twice the memory -- the views are copied.

*/

#pragma once

#include "EvaluatorSecondary.hh"
#include "State.hh"
#include "Factory.hh"

namespace Amanzi {

class EvaluatorSecondaryVectorAsPatch : public EvaluatorSecondary {
 public:
  using EvaluatorSecondary::EvaluatorSecondary;

  virtual void EnsureCompatibility(State& S) override;

  static const std::string eval_type;
  virtual std::string getType() const override { return eval_type; }

 protected:
  // These do the actual work
  virtual void Update_(State& S) override;
  virtual void UpdateDerivative_(State& S, const Key& wrt_key, const Tag& wrt_tag) override {}

 private:
  static Utils::RegisteredFactory<Evaluator, EvaluatorSecondaryVectorAsPatch> reg_;
};

} // namespace Amanzi
