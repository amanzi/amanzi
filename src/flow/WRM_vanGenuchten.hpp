/*
This is the flow component of the Amanzi code. 

Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided in the top-level COPYRIGHT file.

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
  explicit WRM_vanGenuchten(std::string region, double m_, double alpha_, double sr_, double pc0_ = 0.0);
  ~WRM_vanGenuchten() {};
  
  // required methods from the base class
  double k_relative(double pc);
  double saturation(double pc);
  double dSdPc(double pc);  
  double capillaryPressure(double saturation);
  double residualSaturation() { return sr; }

 private:
  double m;  // van Genuchten parameters: m, n, alpha
  double n; 
  const double alpha; 
  const double sr;  // van Genuchten residual saturation

  const double pc0;  // regularization threshold (ususally 0 to 500 Pa)
  double factor_dSdPc;  // frequently used constant
};

}  // namespace AmanziFlow
}  // namespace Amanzi
 
#endif
