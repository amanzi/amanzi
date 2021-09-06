/*
  EOS
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Konstantin Lipnikov (lipnikov@lanl.gov)

  Thermal conductivity for liquid water.
*/

#ifndef AMANZI_EOS_H2O_THERMAL_CONDUCTIVITY_HH_
#define AMANZI_EOS_H2O_THERMAL_CONDUCTIVITY_HH_

#include "Teuchos_ParameterList.hpp"

namespace Amanzi {
namespace AmanziEOS {

// Equation of State model
class H2O_ThermalConductivity {
 public:
  explicit
  H2O_ThermalConductivity(Teuchos::ParameterList& eos_plist);

  virtual double ThermalConductivity(double T);
  virtual double DThermalConductivityDT(double T);

  // error messages  FIXME (move to factory that includes LookupTable)
  int error_code() { return ierr_; }
  std::string error_msg() { return error_msg_; }

 protected:
  virtual void InitializeFromPlist_();

 protected:
  Teuchos::ParameterList eos_plist_;

  // constants for water, hard-coded
  double ka0_, ka1_, ka2_;
  double kref_, Tref_;

  int ierr_;
  std::string error_msg_;
};

}  // namespace AmanziEOS
}  // namespace Amanzi

#endif
