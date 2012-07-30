/*
This is the flow component of the Amanzi code. 
License: BSD
Authors: Neil Carlson (version 1)
         Konstantin Lipnikov (version 2) (lipnikov@lanl.gov)
*/

#ifndef __WATER_RETENTION_MODEL_HPP__
#define __WATER_RETENTION_MODEL_HPP__

#include <string>

namespace Amanzi {
namespace AmanziFlow {

class WaterRetentionModel {
 public:
  virtual double k_relative(double pc) = 0;
  virtual double saturation(double pc) = 0;
  virtual double dSdPc(double pc) = 0;  // derivative of saturation w.r.t. to capillary pressure
  virtual double capillaryPressure(double s) = 0;
  virtual double residualSaturation() = 0;
  
  const std::string region() { return region_; };

 protected:
  void set_region(std::string region) { region_ = region; };
  std::string region_;
};

}  // namespace AmanziFlow
}  // namespace Amanzi
  
#endif
  
