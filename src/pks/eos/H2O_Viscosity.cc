/*
  EOS
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)

  EOS for liquid water. See the permafrost physical properties notes for
  references and documentation of this EOS at:
*/

#include "errors.hh"
#include "H2O_Viscosity.hh"

namespace Amanzi {
namespace AmanziEOS {

H2O_Viscosity::H2O_Viscosity(Teuchos::ParameterList& eos_plist)
  : EOS_Viscosity(eos_plist),
    kav1_(998.333),
    kbv1_(-8.1855),
    kcv1_(0.00585),
    kbv2_(1.3272),
    kcv2_(-0.001053),
    kT1_(293.15) {};


double H2O_Viscosity::Viscosity(double T, double p)
{
  double visc(0.0);

  if (T > 200.0) {
    double A, xi, dT(kT1_ - T);

    if (T < kT1_) {
      A = kav1_ + (kbv1_ + kcv1_ * dT) * dT;
      xi = 1301.0 * (1.0/A - 1.0/kav1_);
    } else {
      A = (kbv2_ + kcv2_ * dT) * dT;
     xi = A / (T - 168.15);
    }
    visc = 0.001 * pow(10.0, xi);
  }

  ierr_ = 0;
  if (visc < 1.0e-16) {
    ierr_ = 1;
    std::stringstream ss;
    ss << "Invalid T=" << T << ", viscosity=" << visc;
    error_msg_ = ss.str();
  }
  return visc;
};


double H2O_Viscosity::DViscosityDT(double T, double p) {
  Errors::Message message("EOS viscosity of water: derivative not implemented");
  Exceptions::amanzi_throw(message);
  return -1.0;
};


double H2O_Viscosity::DViscosityDp(double T, double p) {
  Errors::Message message("EOS viscosity of water: derivative not implemented");
  Exceptions::amanzi_throw(message);
  return -1.0;
};

}  // namespace AmanziEOS
}  // namespace Amanzi
