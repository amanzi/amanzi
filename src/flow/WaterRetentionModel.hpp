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
  virtual double d_saturation(double pc) = 0;
  virtual double capillaryPressure(double sl) = 0;
  
  const std::string region() { return reg; };

 protected:
  void set_region( const std::string region_ ) { reg = region_; };
  std::string reg;
};

}  // namespace AmanziFlow
}  // namespace Amanzi
  
#endif
  
