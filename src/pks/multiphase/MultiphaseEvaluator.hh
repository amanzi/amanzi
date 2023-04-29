/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Multiphase PK

  Adds one function to the secondary evaluator.
*/

#ifndef AMANZI_MULTIPHASE_EVALUATOR_HH_
#define AMANZI_MULTIPHASE_EVALUATOR_HH_

// Amanzi
#include "EvaluatorSecondaryMonotype.hh"
#include "VaporLiquidFactory.hh"

namespace Amanzi {
namespace Multiphase {

class MultiphaseEvaluator
  : public EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace> {
 public:
  MultiphaseEvaluator(Teuchos::ParameterList& plist)
    : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(plist), n_(0){};

  // inteface functions to FieldEvaluator
  MultiphaseEvaluator(const MultiphaseEvaluator& other)
    : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(other)
  {
    n_ = other.n_;
    vapor_liquid_ = other.vapor_liquid_;
  }

  using EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>::Update;
  bool Update(State& S, Key request, bool force)
  {
    bool ok = EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>::Update(S, request);
    if (force) {
      Update_(S);
      requests_.clear();
      requests_.insert(request);
    }
    return ok || force;
  }

  // added interface (WIP)
  // -- modifier
  virtual void set_subvector(int ifield, int n, const std::string& name,
                             Teuchos::ParameterList& plist)
  {
    n_ = n;
    AmanziEOS::VaporLiquidFactory factory(plist);
    vapor_liquid_ = factory.Create(name);
  }

 protected:
  int n_;
  std::shared_ptr<AmanziEOS::VaporLiquid> vapor_liquid_;
};

} // namespace Multiphase
} // namespace Amanzi

#endif
