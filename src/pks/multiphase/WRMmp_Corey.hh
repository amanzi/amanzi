/*
  This is the multiphase component of the Amanzi code. 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Quan Bui (mquanbui@math.umd.edu)

  We use this class to remove effective water retention model.
*/

#ifndef AMANZI_MULTIPHASE_WRM_COREY_HH_
#define AMANZI_MULTIPHASE_WRM_COREY_HH_

#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"

#include "WRMmp.hh"

namespace Amanzi {
namespace Multiphase {

class WRMmp_Corey : public WRMmp {
 public:
  WRMmp_Corey(Teuchos::ParameterList& plist);
  ~WRMmp_Corey(){};

  // required methods from the base class
  virtual double k_relative(double sl, int phase);
  virtual double capillaryPressure(double saturation);
  virtual double dPc_dS(double sl);
  virtual double dKdS(double sl, int phase);

 private:
  double srl_, srg_, pcap_;

  static Utils::RegisteredFactory<WRMmp, WRMmp_Corey> factory_;
};

} // namespace Multiphase
} // namespace Amanzi

#endif
