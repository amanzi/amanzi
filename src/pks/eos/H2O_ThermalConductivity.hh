/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  EOS

  Thermal conductivity for liquid water.
*/

#ifndef AMANZI_EOS_H2O_THERMAL_CONDUCTIVITY_HH_
#define AMANZI_EOS_H2O_THERMAL_CONDUCTIVITY_HH_

#include "Teuchos_ParameterList.hpp"

#include "Factory.hh"

#include "EOS_ThermalConductivity.hh"

namespace Amanzi {
namespace AmanziEOS {

// Equation of State model
class H2O_ThermalConductivity : public EOS_ThermalConductivity {
 public:
  explicit H2O_ThermalConductivity(Teuchos::ParameterList& plist);
  virtual ~H2O_ThermalConductivity(){};

  virtual double ThermalConductivity(double T, double p);
  virtual double DThermalConductivityDT(double T, double p);
  virtual double DThermalConductivityDP(double T, double p) { AMANZI_ASSERT(false); return 0.0; }

  // error messages  FIXME (move to factory that includes LookupTable)
  int error_code() { return ierr_; }
  std::string error_msg() { return error_msg_; }

 protected:
  virtual void InitializeFromPlist_();

 protected:
  // constants for water, hard-coded
  double ka0_, ka1_, ka2_;
  double kref_, Tref_;

  int ierr_;
  std::string error_msg_;

 private:
  static Utils::RegisteredFactory<EOS_ThermalConductivity, H2O_ThermalConductivity> reg_;
};

} // namespace AmanziEOS
} // namespace Amanzi

#endif
