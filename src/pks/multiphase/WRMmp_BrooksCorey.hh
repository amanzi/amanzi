/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Quan Bui (mquanbui@math.umd.edu)
*/

/*
  This is the multiphase component of the Amanzi code.

  We use this class to set up simple water retention model
  for decoupled multiphase flow
*/

#ifndef AMANZI_MULTIPHASE_WRM_BROOKS_COREY_HH_
#define AMANZI_MULTIPHASE_WRM_BROOKS_COREY_HH_

#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"

#include "WRMmp.hh"

namespace Amanzi {
namespace Multiphase {

class WRMmp_BrooksCorey : public WRMmp {
 public:
  WRMmp_BrooksCorey(Teuchos::ParameterList& plist);
  ~WRMmp_BrooksCorey(){};

  // required methods from the base class
  virtual double k_relative(double Sw, int phase);
  virtual double capillaryPressure(double saturation);
  virtual double dPc_dS(double saturation);
  virtual double dKdS(double Sw, int phase);

 private:
  void Init_(double S_rw, double S_rn, double pd, double lambda);

 private:
  double S_rw_, S_rn_, pd_, lambda_;

  static Utils::RegisteredFactory<WRMmp, WRMmp_BrooksCorey> reg_;
};

} // namespace Multiphase
} // namespace Amanzi

#endif
