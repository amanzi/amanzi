/*
This is the flow component of the Amanzi code. 
License: BSD
Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef __VAN_ANALYTIC_MODEL_HPP__
#define __VAN_ANALYTIC_MODEL_HPP__

#include "WaterRetentionModel.hpp"

/*
 We use this class to test convergence of discretization schemes.
 It uses the simplest model for relrative permeability, k_rel = 1 / (1 + p^2).
*/

namespace Amanzi {
namespace AmanziFlow {

class WRM_analytic : public WaterRetentionModel {
 public:
  WRM_analytic(std::string region_);
  ~WRM_analytic() {};
  
  // requires methods from the base class
  double k_relative(double p);
  double saturation(double p);
  double d_saturation(double p);  
  double pressure(double saturation);

 private:
  double m, n, alpha, atm_pressure;
};

}  // namespace AmanziFlow
}  // namespace Amanzi
 
#endif
