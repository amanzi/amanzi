/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  ATS

  ViscosityModel for liquid water.  See the permafrost physical properties notes for
  references and documentation of this EOS at:

  http://software.lanl.gov/ats/trac

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef FLOWRELATIONS_EOS_WATER_HH_
#define FLOWRELATIONS_EOS_WATER_HH_

#include "Teuchos_ParameterList.hpp"

#include "factory.hh"
#include "viscosity_model.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

// Equation of State model
class ViscosityModelWater : public ViscosityModel {

public:
  explicit
  ViscosityModelWater(Teuchos::ParameterList& eos_plist);

  ViscosityModelWater(const ViscosityModelWater& other);

  virtual Teuchos::RCP<FieldModel> Clone() const;

  virtual double Viscosity(double T);
  virtual double DViscosityDT(double T);

private:
  virtual void InitializeFromPlist_();

  // constants for water, hard-coded because it would be crazy to try to come
  // up with names for these in a parameter list...
  // -- temperature dependence of viscosity < T1
  const double kav1_, kbv1_, kcv1_;

  // -- temperature dependence of viscosity > T1
  const double kbv2_, kcv2_, kT1_;

  static Utils::RegisteredFactory<FieldModel,ViscosityModelWater> factory_;

};

} // namespace
} // namespace
} // namespace

#endif
