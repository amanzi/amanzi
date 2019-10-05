/*
  This is the multiphase component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Quan Bui (mquanbui@math.umd.edu)

  We use this class to set up simple water retention model 
  for decoupled multiphase flow
*/

#ifndef AMANZI_SIMPLE_MODEL_HH_
#define AMANZI_SIMPLE_MODEL_HH_

#include "WaterRetentionModel.hh"

namespace Amanzi {
namespace Multiphase {

class WRM_Simple : public WaterRetentionModel {
 public:
  explicit WRM_Simple(std::string region, double S_rw, double S_rn, double coef);
  ~WRM_Simple() {};
  
  // required methods from the base class
  double k_relative(double Sw, std::string phase_name);
  double capillaryPressure(double saturation);
  double dPc_dS(double saturation);
  double residualSaturation(std::string phase_name);
  double dKdS(double Sw, std::string phase_name);

 private:
  double S_rw_, S_rn_, coef_, exponent_;
};

}  // namespace Flow
}  // namespace Amanzi
 
#endif
