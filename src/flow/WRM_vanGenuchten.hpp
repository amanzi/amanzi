/*
This is the flow component of the Amanzi code. 
License: BSD
Authors: Neil Carlson (version 1)
         Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
*/

#ifndef __VAN_GENUCHTEN_MODEL_HPP__
#define __VAN_GENUCHTEN_MODEL_HPP__

#include "WaterRetentionModel.hpp"

namespace Amanzi {
namespace AmanziFlow {

class WRM_vanGenuchten : public WaterRetentionModel {
 public:
  WRM_vanGenuchten(std::string region_, double m_, double alpha_, double sr_);
  ~WRM_vanGenuchten() {};
  
  // requires methods from the base class
  double k_relative(double pc);
  double saturation(double pc);
  double d_saturation(double pc);  
  double capillaryPressure(double saturation);

 private:
  double m;  // van Genuchten parameters: m, n, alpha
  double n; 
  const double alpha; 
  const double sr;  // van Genuchten effective saturation
};

}  // namespace AmanziFlow
}  // namespace Amanzi
 
#endif
