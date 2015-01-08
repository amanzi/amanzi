/*
  This is the EOS component of the ATS and Amanzi codes.
   
  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)

  Viscosity for liquid water.  See the permafrost physical properties notes for
  references and documentation of this EOS at:
*/

#ifndef AMANZI_EOS_VISCOSITY_WATER_HH_
#define AMANZI_EOS_VISCOSITY_WATER_HH_

#include "Teuchos_ParameterList.hpp"

#include "factory.hh"
#include "dbc.hh"
#include "viscosity_relation.hh"

namespace Amanzi {
namespace Relations {

// Equation of State model
class ViscosityWater : public ViscosityRelation {
 public:
  explicit
  ViscosityWater(Teuchos::ParameterList& eos_plist);

  virtual double Viscosity(double T);
  virtual double DViscosityDT(double T);

 protected:
  Teuchos::ParameterList eos_plist_;

  // constants for water, hard-coded because it would be crazy to try to come
  // up with names for these in a parameter list...
  // -- temperature dependence of viscosity < T1
  const double kav1_, kbv1_, kcv1_;

  // -- temperature dependence of viscosity > T1
  const double kbv2_, kcv2_, kT1_;

 private:
  static Utils::RegisteredFactory<ViscosityRelation,ViscosityWater> factory_;
};

}  // namespace Relations
}  // namespace Amanzi

#endif
