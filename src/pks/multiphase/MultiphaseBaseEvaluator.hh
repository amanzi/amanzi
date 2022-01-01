/*
  Multiphase PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)

  Adds one function to the second variable evaluator.
*/

#ifndef AMANZI_MULTIPHASE_BASE_EVALUATOR_HH_
#define AMANZI_MULTIPHASE_BASE_EVALUATOR_HH_

// Amanzi
#include "EvaluatorSecondaryMonotype.hh"

namespace Amanzi {
namespace Multiphase {

class MultiphaseBaseEvaluator : public EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace> {
 public:
  MultiphaseBaseEvaluator(Teuchos::ParameterList& plist)
    : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(plist),
      n_(0), kH_(1.0) {};

  // inteface functions to FieldEvaluator
  MultiphaseBaseEvaluator(const MultiphaseBaseEvaluator& other)
      : EvaluatorSecondaryMonotype<CompositeVector, CompositeVectorSpace>(other) {
    n_ = other.n_; kH_ = other.kH_;
  }

  bool Update(State& S, Key request, bool force) {
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
  virtual void set_subvector(int ifield, int n, double kH) { n_ = n; kH_ = kH; } 

 protected:
  int n_;
  double kH_;
};

}  // namespace Multiphase
}  // namespace Amanzi

#endif
