/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Mechanics PK

  Evaluator for computing hydrostatic stress.
*/

#include "PDE_Elasticity.hh"

#include "ShearStrainEvaluator.hh"

namespace Amanzi {
namespace Mechanics {

/* ******************************************************************
* Constructor
****************************************************************** */
ShearStrainEvaluator::ShearStrainEvaluator(Teuchos::ParameterList& plist)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(plist)
{
  displacement_key_ = plist.get<std::string>("displacement key", "displacement");
  dependencies_.insert(std::make_pair(displacement_key_, Tags::DEFAULT));
}


/* ******************************************************************
* A copy constructor.
****************************************************************** */
ShearStrainEvaluator::ShearStrainEvaluator(const ShearStrainEvaluator& other)
  : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(other),
    displacement_key_(other.displacement_key_),
    op_(other.op_)
{}


/* ******************************************************************
* Clone with unclear yet purpose.
****************************************************************** */
Teuchos::RCP<Evaluator>
ShearStrainEvaluator::Clone() const
{
  return Teuchos::rcp(new ShearStrainEvaluator(*this));
}


/* ******************************************************************
* Evaluator
****************************************************************** */
void
ShearStrainEvaluator::Evaluate_(const State& S, const std::vector<CompositeVector*>& results)
{
  if (op_.get()) {
    auto u = S.Get<CompositeVector>(displacement_key_, Tags::DEFAULT);
    u.ScatterMasterToGhosted();

    auto& shear_c = *results[0]->ViewComponent("cell");
    int ncells = shear_c.MyLength();


    for (int c = 0; c < ncells; ++c) {
      auto Tc = op_->ComputeCellStrain(u, c);

      int d = Tc.dimension();
      double trace = Tc.Trace() / d;
      for (int i = 0; i < d; ++i) Tc(i, i) -= trace;

      double gamma(0.0);
      for (int i = 0; i < d; ++i) {
        gamma += Tc(i, i) * Tc(i, i);
        for (int j = i + 1; j < d; ++j) gamma += 2 * Tc(i, j) * Tc(i, j);
      }
      shear_c[0][c] = std::sqrt(2 * gamma);
    }
  } else {
    // operator may not be available during initialization time
    results[0]->PutScalar(0.0);
  }
}

} // namespace Mechanics
} // namespace Amanzi
