/*
  EOS
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)

  Constant prescribed viscosity.
*/

#ifndef AMANZI_EOS_H2O_VISCOSITY_CONSTANT_HH_
#define AMANZI_EOS_H2O_VISCOSITY_CONSTANT_HH_

#include "Teuchos_ParameterList.hpp"

#include "dbc.hh"
#include "Factory.hh"

#include "EOS_Viscosity.hh"

namespace Amanzi {
namespace AmanziEOS {

class H2O_ViscosityConstant : public EOS_Viscosity {
 public:
  explicit H2O_ViscosityConstant(Teuchos::ParameterList& visc_plist);

  virtual double Viscosity(double T, double p) override { return visc_; }
  virtual double DViscosityDT(double T, double p) override { return 0.0; }
  virtual double DViscosityDp(double T, double p) override { return 0.0; }

 protected:
  virtual void InitializeFromPlist_();

  double visc_;

 private:
  static Utils::RegisteredFactory<EOS_Viscosity, H2O_ViscosityConstant> factory_;
};

} // namespace AmanziEOS
} // namespace Amanzi

#endif
