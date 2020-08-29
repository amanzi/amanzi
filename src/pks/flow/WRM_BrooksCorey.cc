/*
  Flow PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

#include <cmath>
#include <string>

#include "errors.hh"
#include "FlowDefs.hh"
#include "WRM_BrooksCorey.hh"

namespace Amanzi {
namespace Flow {

/* ******************************************************************
* Setup fundamental parameters for this model.
* Default value of the regularization interval is pc0 = 0.                                           
****************************************************************** */
WRM_BrooksCorey::WRM_BrooksCorey(Teuchos::ParameterList& plist)
{
  double lambda = plist.get<double>("Brooks Corey lambda", FLOW_WRM_EXCEPTION);
  double alpha = plist.get<double>("Brooks Corey alpha", FLOW_WRM_EXCEPTION);
  double l = plist.get<double>("Brooks Corey l", FLOW_WRM_BROOKS_COREY_L);
  double sr = plist.get<double>("residual saturation liquid", FLOW_WRM_EXCEPTION);
  double pc0 = plist.get<double>("regularization interval", FLOW_WRM_REGULARIZATION_INTERVAL);
  std::string krel_function = plist.get<std::string>("relative permeability model", "Mualem");

  Init_(lambda, l, alpha, sr, krel_function, pc0);
}


/* ******************************************************************
* Setup fundamental parameters for this model.
* Default value of the regularization interval is pc0 = 0.                                           
****************************************************************** */
void WRM_BrooksCorey::Init_(
    double lambda, double l, double alpha, 
    double sr, std::string& krel_function, double pc0)
{
  lambda_ = lambda;
  l_ = l; 
  alpha_ = alpha;
  sr_ = sr;
  pc0_ = pc0;

  Errors::Message msg;
  if (l_ < 0.0 || lambda_ < 0.0 || sr_ < 0.0 || pc0_ < 0.0) {
    msg << "Brooks Corey: negative parameter in a water retention model.";
    Exceptions::amanzi_throw(msg);
  }
  if (sr_ > 1.0) {
    msg << "Brooks Corey: residual saturation is greater than 1.";
    Exceptions::amanzi_throw(msg);
  }

  if (krel_function == "Mualem") {
    factor_ = -2.0 - (l_ + 2.0) * lambda_;
  } else if (krel_function == "Burdine") {
    factor_ = -2.0 - 3.0 * lambda_;
  }

  pc_bubble_ = 1.0 / alpha_;
  a_ = b_ = 0;
  
  if (pc0_ > 0.0) {
    pc0_ += pc_bubble_;
    double k0 = k_relative(pc0_) - 1.0;
    double k0p = dKdPc(pc0_);
    double pc0_2 = pc0 * pc0;
    double pc0_3 = pc0_2 * pc0;

    a_ = (3 * k0 - k0p * pc0) / pc0_2;
    b_ = (k0p * pc0 - 2 * k0) / pc0_3;
  }
}


/* ******************************************************************
* Relative permeability formula: input is capillary pressure pc.
* The original curve is regulized on interval (0, pc0) using the 
* Hermite interpolant of order 3. Formulas (3.14)-(3.15).   
****************************************************************** */
double WRM_BrooksCorey::k_relative(double pc) const
{
  if (pc <= pc_bubble_) {
    return 1.0;
  } else if (pc >= pc0_) {
    return pow(alpha_ * pc, factor_);
  } else {
    double pc_2 = (pc - pc_bubble_) * (pc - pc_bubble_);
    double pc_3 = pc_2 * (pc - pc_bubble_);
    return 1.0 + a_ * pc_2 + b_ * pc_3;
  }
}


/* ******************************************************************
* Saturation formula (3.5)-(3.8).                                         
****************************************************************** */
double WRM_BrooksCorey::saturation(double pc) const
{
  if (pc > pc_bubble_) {
    return pow(alpha_ * pc, -lambda_) * (1.0 - sr_) + sr_;
  } else {
    return 1.0;
  }
}


/* ******************************************************************
* Derivative of the saturation formula w.r.t. capillary pressure.
* Warning: remember that dSdP = -dSdPc.                                        
****************************************************************** */
double WRM_BrooksCorey::dSdPc(double pc) const
{
  if (pc > pc_bubble_) {
    return -pow(alpha_ * pc, -lambda_ - 1.0) * (1.0 - sr_) * alpha_ * lambda_;
  } else {
    return 0.0;
  }
}


/* ******************************************************************
* Pressure as a function of saturation, formula (3.9).                                       
****************************************************************** */
double WRM_BrooksCorey::capillaryPressure(double s) const
{
  if (s == sr_) return 0.0;
  double se = (s - sr_) / (1.0 - sr_);
  return pow(se, -1.0/lambda_) / alpha_;
}


/* ******************************************************************
* Derivative of the original relative permeability w.r.t. capillary pressure.                                     
****************************************************************** */
double WRM_BrooksCorey::dKdPc(double pc)  const
{
  if (pc > pc_bubble_) {
    return factor_ * alpha_ * pow(alpha_ * pc, factor_ - 1.0);
  } else {
    return 0.0;
  }
}

}  // namespace Flow
}  // namespace Amanzi

