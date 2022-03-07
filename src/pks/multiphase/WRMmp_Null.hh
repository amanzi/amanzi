/*
  This is the multiphase component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Quan Bui (mquanbui@math.umd.edu)

  We use this class to remove effective water retention model.
*/

#ifndef AMANZI_MULTIPHASE_WRM_NULL_HH_
#define AMANZI_MULTIPHASE_WRM_NULL_HH_

#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"

#include "WRMmp.hh"

namespace Amanzi {
namespace Multiphase {

class WRMmp_Null : public WRMmp {
 public:
  WRMmp_Null(Teuchos::ParameterList& plist) {};
  ~WRMmp_Null() {};
  
  // required methods from the base class
  virtual double k_relative(double Sw, int phase) { return 1.0; }
  virtual double capillaryPressure(double saturation) { return 0.0; }
  virtual double dPc_dS(double saturation) { return 0.0; }
  virtual double dKdS(double Sw, int phase) { return 0.0; }

 private:
  static Utils::RegisteredFactory<WRMmp, WRMmp_Null> factory_;
};

}  // namespace Multiphase
}  // namespace Amanzi
 
#endif
