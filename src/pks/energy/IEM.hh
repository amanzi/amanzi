/*
  Energy

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon

  Internal energy model is function of temperature only.
  UNITS: J/{mol/kg}
*/

#ifndef AMANZI_ENERGY_IEM_HH_
#define AMANZI_ENERGY_IEM_HH_

// #include "Teuchos_ParameterList.hpp"

namespace Amanzi {
namespace Energy {

class IEM {
 public:
  IEM() : ierr_(0){};
  virtual ~IEM() {}

  // IEM(Teuchos::ParameterList& plist);
  virtual double InternalEnergy(double T, double p) = 0;
  virtual double DInternalEnergyDT(double T, double p) = 0;
  virtual double DInternalEnergyDp(double T, double p) = 0;

  // error messages
  int error_code() { return ierr_; }
  std::string error_msg() { return error_msg_; }

 protected:
  int ierr_;
  std::string error_msg_;
};

} // namespace Energy
} // namespace Amanzi

#endif
