/*
  EOS
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Author: Ethan Coon (ecoon@lanl.gov)

  Viscosity for liquid water.
*/

#ifndef AMANZI_EOS_VISCOSITY_WATER_HH_
#define AMANZI_EOS_VISCOSITY_WATER_HH_

#include "Teuchos_ParameterList.hpp"

#include "dbc.hh"
#include "Factory.hh"

#include "Viscosity_Base.hh"

namespace Amanzi {
namespace AmanziEOS {

// Equation of State model
class Viscosity_Water : public Viscosity_Base {
 public:
  explicit
  Viscosity_Water(Teuchos::ParameterList& eos_plist);

  virtual double Viscosity(double T);
  virtual double DViscosityDT(double T);

 protected:
  Teuchos::ParameterList eos_plist_;

  // constants for water, hard-coded 
  const double kav1_, kbv1_, kcv1_;
  const double kbv2_, kcv2_, kT1_;

 private:
  static Utils::RegisteredFactory<Viscosity_Base, Viscosity_Water> factory_;
};

}  // namespace AmanziEOS
}  // namespace Amanzi

#endif
