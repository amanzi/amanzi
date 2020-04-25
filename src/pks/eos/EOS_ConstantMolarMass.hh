/*
  EOS
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)

  EOS_ConstantMolarMass -- intermediate class for default implementations
  of EOS with constant molar masses.

  Note, while instantiating this class does work, it should not be done.
  Instead, this class is intended to be inherited and either the Molar or Mass
  methods replaced.  Use as it stands results in an infinite recursion...
*/

#ifndef AMANZI_EOS_CONSTANT_MOLAR_MASS_HH_
#define AMANZI_EOS_CONSTANT_MOLAR_MASS_HH_

#include "EOS.hh"

namespace Amanzi {
namespace AmanziEOS {

class EOS_ConstantMolarMass : public EOS {
 public:
  EOS_ConstantMolarMass() : M_(0.0) {};
  explicit EOS_ConstantMolarMass(double M) : M_(M) {};

  virtual double MolarDensity(double T, double p) {
    return MassDensity(T, p) / M_;
  }

  virtual double DMolarDensityDT(double T, double p) {
    return DMassDensityDT(T, p) / M_;
  }

  virtual double DMolarDensityDp(double T, double p) {
    return DMassDensityDp(T, p) / M_;
  }

  virtual double MassDensity(double T, double p) {
    return MolarDensity(T, p) * M_;
  }

  virtual double DMassDensityDT(double T, double p) {
    return DMolarDensityDT(T, p) * M_;
  }

  virtual double DMassDensityDp(double T, double p) {
    return DMolarDensityDp(T, p) * M_;
  }

  virtual bool IsConstantMolarMass() { return true; }
  virtual double MolarMass() { return M_; }

 protected:
  double M_;
};

}  // namespace AmanziEOS
}  // namespace Amanzi

#endif
