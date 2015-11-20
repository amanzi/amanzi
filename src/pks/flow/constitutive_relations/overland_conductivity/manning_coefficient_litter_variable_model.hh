/*
  The manning coefficient with variable litter model is an algebraic model with dependencies.
  
  Constant values.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef AMANZI_FLOW_MANNING_LITTER_COEFFICIENT_VARIABLE_MODEL_HH_
#define AMANZI_FLOW_MANNING_LITTER_COEFFICIENT_VARIABLE_MODEL_HH_

#include "manning_coefficient_litter_model.hh"

namespace Amanzi {
namespace Flow {
namespace Relations {

class ManningCoefficientLitterVariableModel : public ManningCoefficientLitterModel {
 public:
  ManningCoefficientLitterVariableModel(Teuchos::ParameterList& plist) {
    n_bg_ = plist.get<double>("manning coefficient bare ground [s * m^-1/3]", 0.02);
    n_l_ = plist.get<double>("manning coefficient litter [s * m^-1/3]", 0.1);
  }
  
  double ManningCoefficient(double ld, double pd) const {
    return (pd <= ld) ? n_l_ : (ld*n_l_ + n_bg_*(-ld + pd))/pd;
  }
    

  double DManningCoefficientDLitterThickness(double ld, double pd) const {
    return (pd <= ld) ? 0 : (-n_bg_ + n_l_)/pd;
  }


  double DManningCoefficientDPondedDepth(double ld, double pd) const {
    return (pd <= ld) ? 0 : n_bg_/pd - (ld*n_l_ + n_bg_*(-ld + pd))/std::pow(pd, 2);
  }

 protected:
  double n_bg_;
  double n_l_;
  
};

} //namespace
} //namespace
} //namespace

#endif
