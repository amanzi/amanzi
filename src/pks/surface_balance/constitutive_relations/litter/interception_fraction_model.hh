/*
  The interception fraction model is an algebraic model with dependencies.

  Generated via evaluator_generator with:

    
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_SURFACEBALANCE_INTERCEPTION_FRACTION_MODEL_HH_
#define AMANZI_SURFACEBALANCE_INTERCEPTION_FRACTION_MODEL_HH_

namespace Amanzi {
namespace SurfaceBalance {
namespace Relations {

class InterceptionFractionModel {

 public:
  explicit
  InterceptionFractionModel(Teuchos::ParameterList& plist);

  double InterceptionFraction(double ai) const;

  double DInterceptionFractionDAreaIndex(double ai) const;
  
 protected:
  void InitializeFromPlist_(Teuchos::ParameterList& plist);

 protected:

  double alpha_;

};

} //namespace
} //namespace
} //namespace

#endif
