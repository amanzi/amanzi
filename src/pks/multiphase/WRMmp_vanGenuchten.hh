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

#ifndef AMANZI_MULTIPHASE_WRM_VAN_GENUCHTEN_HH_
#define AMANZI_MULTIPHASE_WRM_VAN_GENUCHTEN_HH_

#include "Teuchos_ParameterList.hpp"

// Amanzi
#include "Factory.hh"
#include "SplinePolynomial.hh"

// Multiphase
#include "WRMmp.hh"

namespace Amanzi {
namespace Multiphase {

class WRMmp_vanGenuchten : public WRMmp {
 public:
  WRMmp_vanGenuchten(Teuchos::ParameterList& plist);
  ~WRMmp_vanGenuchten() {};
  
  // required methods from the base class
  virtual double k_relative(double sl, int phase);
  virtual double capillaryPressure(double saturation);
  virtual double dPc_dS(double saturation);
  virtual double dKdS(double sl, int phase);

 private:
  void Init_(double srw, double srn, double n, double Pr, double reg_kr, double reg_pc);

  double k_relative_gas_(double sle);
  double k_relative_liquid_(double sle);

  double dKdS_gas_(double sle);
  double dKdS_liquid_(double sle);

  double capillaryPressure_(double sle);
  double dPc_dS_(double sle);

 private:
  double Pr_, srl_, srg_, n_, m_;

  double reg_kl_;
  WhetStone::SplineCubic spline_kl_;

  double reg_pc_;
  WhetStone::SplineExteriorLinear spline_pc_;

  static Utils::RegisteredFactory<WRMmp, WRMmp_vanGenuchten> factory_;
};

}  // namespace Multiphase
}  // namespace Amanzi
 
#endif
