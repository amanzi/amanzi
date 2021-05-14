/*
  EOS
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)

  EOS is purely virtual base class for an equation of state.
*/

#ifndef AMANZI_EOS_HH_
#define AMANZI_EOS_HH_

#include "Teuchos_ParameterList.hpp"

namespace Amanzi {
namespace AmanziEOS {

class EOS_Density {
 public:
  EOS_Density(Teuchos::ParameterList& eos_plist)
    : eos_plist_(eos_plist) {
    InitializeFromPlist_();
  }
  virtual ~EOS_Density() {};

  // Virtual methods that form the EOS
  virtual double Density(double T, double p) = 0;
  virtual double DDensityDT(double T, double p) = 0;
  virtual double DDensityDp(double T, double p) = 0;

  virtual double MolarDensity(double T, double p) = 0;
  virtual double DMolarDensityDT(double T, double p) = 0;
  virtual double DMolarDensityDp(double T, double p) = 0;

  double MolarMass() { return M_; }

 protected:
  virtual void InitializeFromPlist_() {
  M_ = eos_plist_.get<double>("molar mass", 18.0153e-03);

  if (eos_plist_.isParameter("molar density"))
    rho_ = eos_plist_.get<double>("molar density") * M_;
  else 
    rho_ = eos_plist_.get<double>("density", 997.07);
  }

 protected:
  Teuchos::ParameterList eos_plist_;
  double M_, rho_;
};

}  // namespace AmanziEOS
}  // namespace Amanzi

#endif
