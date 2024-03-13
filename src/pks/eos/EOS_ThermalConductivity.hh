/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  EOS

  Equation of state for thermal conductivity = f(T)
*/

#ifndef AMANZI_EOS_THERMAL_CONDUCTIVITY_HH_
#define AMANZI_EOS_THERMAL_CONDUCTIVITY_HH_

#include "Teuchos_ParameterList.hpp"

namespace Amanzi {
namespace AmanziEOS {

class EOS_ThermalConductivity {
 public:
  EOS_ThermalConductivity(Teuchos::ParameterList& plist) : plist_(plist), ierr_(0){};
  virtual ~EOS_ThermalConductivity(){};

  virtual double ThermalConductivity(double T, double p) = 0;
  virtual double DThermalConductivityDT(double T, double p) = 0;
  virtual double DThermalConductivityDP(double T, double p) = 0;

  // error messages
  int error_code() { return ierr_; }
  std::string error_msg() { return error_msg_; }

 protected:
  Teuchos::ParameterList plist_;

  int ierr_;
  std::string error_msg_;
};

} // namespace AmanziEOS
} // namespace Amanzi

#endif
