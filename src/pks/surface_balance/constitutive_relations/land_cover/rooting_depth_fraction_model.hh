/*
  The rooting depth fraction model is an algebraic model with dependencies.

  Generated via evaluator_generator with:
Rooting depth function.

Sets the root fraction as a function of depth,

F_root =  ( a*exp(-az) + b*exp(-bz) ) / 2

This function is such that the integral over depth = [0,inf) is 1, but
an artificial cutoff is generated.


    
  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_FLOW_ROOTING_DEPTH_FRACTION_MODEL_HH_
#define AMANZI_FLOW_ROOTING_DEPTH_FRACTION_MODEL_HH_

namespace Amanzi {
namespace LandCover {
namespace Relations {

class RootingDepthFractionModel {

 public:
  explicit
  RootingDepthFractionModel(Teuchos::ParameterList& plist);

  double RootingDepthFraction(double z) const;

  double DRootingDepthFractionDDepth(double z) const;
  
 protected:
  void InitializeFromPlist_(Teuchos::ParameterList& plist);

 protected:

  double a_;
  double b_;
  double z_max_;

};

} //namespace
} //namespace
} //namespace

#endif