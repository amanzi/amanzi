/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.
*/

/*
  This is the multiphase component of the Amanzi code.

  We use this class to set up a custom water retention model used in 1D manufactured solution.
*/

#ifndef AMANZI_MULTIPHASE_WRM_CUSTOM_HH_
#define AMANZI_MULTIPHASE_WRM_CUSTOM_HH_

#include "Factory.hh"

#include "WRMmp.hh"

namespace Amanzi {
namespace Multiphase {

class WRMmp_Custom : public WRMmp {
 public:
  WRMmp_Custom(Teuchos::ParameterList& plist) { double coef = plist.get<double>("coefficient", MULTIPHASE_WRM_EXCEPTION); }
  ~WRMmp_Custom(){};

  // required methods from the base class
  virtual double k_relative(double Sw, int phase) { return 1.0; }
  virtual double capillaryPressure(double saturation) { return 1.0; }
  virtual double dPc_dS(double saturation) { return 0.0; }
  virtual double dKdS(double Sw, int phase) { return 0.0; }

 private:
  void Init_(double coef) { coef = 1.0; }

 private:
  double S_rw_, S_rn_, coef_, exponent_;

  static Utils::RegisteredFactory<WRMmp, WRMmp_Custom> reg_;
};

} // namespace Multiphase
} // namespace Amanzi

#endif
