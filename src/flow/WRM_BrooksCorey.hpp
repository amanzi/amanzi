/*
This is the flow component of the Amanzi code. 

Copyright 2010-2012 held jointly by LANS/LANL, LBNL, and PNNL. 
Amanzi is released under the three-clause BSD License. 
The terms of use and "as is" disclaimer for this license are 
provided in the top-level COPYRIGHT file.

Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#ifndef __BROOKS_COREY_MODEL_HPP__
#define __BROOKS_COREY_MODEL_HPP__

#include "WaterRetentionModel.hpp"

namespace Amanzi {
namespace AmanziFlow {

class WRM_BrooksCorey : public WaterRetentionModel {
 public:
  explicit WRM_BrooksCorey(std::string region, double lambda, double alpha, double sr, 
                           std::string krel_function, double pc0 = 0.0);
  ~WRM_BrooksCorey() {};
  
  // required methods from the base class
  double k_relative(double pc);
  double saturation(double pc);
  double dSdPc(double pc);  
  double capillaryPressure(double saturation);
  double residualSaturation() { return sr_; }

 private:
  double lambda_, alpha_;  // Brooks and Corey parameters: lambda, alpha
  double sr_;  // residual saturation
  int krel_function_;  // Mualem or Burdine

  double pc0_;  // regularization threshold (ususally 0 to 500 Pa)
  double factor_;  // frequently used constant
};

}  // namespace AmanziFlow
}  // namespace Amanzi
 
#endif
