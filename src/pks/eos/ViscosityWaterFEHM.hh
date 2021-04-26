/*
  EOS
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)

  Viscosity for liquid water for T between 0.001 and 360 C from FEHM manual
*/

#ifndef AMANZI_EOS_VISCOSITY_WATER_FEHM_HH_
#define AMANZI_EOS_VISCOSITY_WATER_FEHM_HH_

#include "Teuchos_ParameterList.hpp"

#include "dbc.hh"
#include "Factory.hh"

#include "ViscosityBase.hh"

namespace Amanzi {
namespace AmanziEOS {

// Equation of State model
class ViscosityWaterFEHM : public ViscosityBase {
 public:
  explicit
  ViscosityWaterFEHM(Teuchos::ParameterList& eos_plist);

  virtual double Viscosity(double T, double p);
  virtual double DViscosityDT(double T, double p);
  virtual double DViscosityDp(double T, double p);

 protected:
  double y0_, y1_, y2_, y3_, y4_, y5_, y6_, y7_, y8_, y9_;
  double z0_, z1_, z2_, z3_, z4_, z5_, z6_, z7_, z8_, z9_;
  double T0_;

 private:
  static Utils::RegisteredFactory<ViscosityBase, ViscosityWaterFEHM> factory_;
};

}  // namespace AmanziEOS
}  // namespace Amanzi

#endif
