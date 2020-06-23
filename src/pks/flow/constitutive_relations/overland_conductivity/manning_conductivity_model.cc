/* -*-  mode: c++; indent-tabs-mode: nil -*- */

/*
  Evaluates the conductivity of surface flow as a function of ponded
  depth and surface slope using Manning's model.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#include "manning_conductivity_model.hh"

namespace Amanzi {
namespace Flow {

ManningConductivityModel::ManningConductivityModel(Teuchos::ParameterList& plist) :
    plist_(plist) {

  slope_regularization_ = plist_.get<double>("slope regularization epsilon", 1.e-8);
  manning_exp_ = plist_.get<double>("Manning exponent");
}

double ManningConductivityModel::Conductivity(double depth, double slope, double coef) {
  if (depth <= 0.) return 0.;

  double exponent = manning_exp_ + 1.0;
  double scaling = coef * std::sqrt(std::max(slope, slope_regularization_));
  return std::pow(depth, exponent) / scaling;
}

double ManningConductivityModel::Conductivity(double depth, double slope, double coef, double pdd, double frac_cond, double beta) {
  if (pdd <= 0.) return 0.;

  double exponent = manning_exp_; // 2/3.
  double scaling = coef * std::sqrt(std::max(slope, slope_regularization_));
  double val = std::pow(frac_cond,beta);

  return pdd*frac_cond *val*std::pow(std::max(pdd,0.), exponent) / scaling;

}

double ManningConductivityModel::DConductivityDDepth(double depth, double slope, double coef) {
  if (depth <= 0.) return 0.;
  double exponent = manning_exp_ + 1.0;
  double scaling = coef * std::sqrt(std::max(slope, slope_regularization_));
  return std::pow(std::max(depth,0.), exponent - 1.) * exponent / scaling;
}

double ManningConductivityModel::DConductivityDDepth(double depth, double slope, double coef, double pd_depth, double frac_cond, double beta) {  
  if (pd_depth <= 0.) return 0.;
  //Errors::Message message("Manning Conductivity Model: Derivaritve not implemented for the Subgrid Model."); 
  // Exceptions::amanzi_throw(message);
  std::cout<<"Manning Conductivity Model: Derivaritve not implemented for the Subgrid Model.\n"; abort();  
  double exponent = manning_exp_;
  double scaling = coef * std::sqrt(std::max(slope, slope_regularization_));
  double Dpd1 = depth*std::pow(std::max(pd_depth,0.), exponent - 1.)*exponent;
  double Dpd2 = std::pow(std::max(pd_depth,0.), exponent);
  return 0;
  //return  (Dpd1  + Dpd2) / (scaling * hdmax);
  }




} // namespace
} // namespace
