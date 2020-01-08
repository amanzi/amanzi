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

#include "Factory.hh"

#include "WRMmp.hh"

namespace Amanzi {
namespace Multiphase {

class WRMmp_vanGenuchten : public WRMmp {
 public:
  WRMmp_vanGenuchten(Teuchos::ParameterList& plist);
  ~WRMmp_vanGenuchten() {};
  
  // required methods from the base class
  double k_relative(double sl, std::string phase_name);
  double capillaryPressure(double saturation);
  double dPc_dS(double saturation);
  double dKdS(double sl, std::string phase_name);

  void Init_(double srw, double srn, double n, double Pr);

 private:
  double VGM(double sg);
  double mod_VGM(double sg);
  double deriv_VGM(double sg);
  double deriv_mod_VGM(double sg);

  double Pr_, srl_, srg_, n_, m_, eps_;

  static Utils::RegisteredFactory<WRMmp, WRMmp_vanGenuchten> factory_;
};

}  // namespace Multiphase
}  // namespace Amanzi
 
#endif
