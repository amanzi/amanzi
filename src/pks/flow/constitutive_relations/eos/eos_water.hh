/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  ATS

  EOS for liquid water.  See the permafrost physical properties notes for
  references and documentation of this EOS at:

  http://software.lanl.gov/ats/trac

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef FLOWRELATIONS_EOS_WATER_HH_
#define FLOWRELATIONS_EOS_WATER_HH_

#include "Teuchos_ParameterList.hpp"

#include "factory_with_state.hh"
#include "eos.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

// Equation of State model
class EOSWater : public EOS {

public:
  EOSWater(Teuchos::ParameterList& eos_plist, const Teuchos::Ptr<State>& S);
  EOSWater(const EOSWater& other);

  virtual Teuchos::RCP<FieldModel> Clone() const;

  virtual double Density(double T, double p);
  virtual double DDensityDT(double T, double p);
  virtual double DDensityDp(double T, double p);

  virtual double molar_mass() { return M_; }
  virtual bool is_molar_basis() { return false; }

private:
  virtual void InitializeFromPlist_();

  double M_;

  // constants for water, hard-coded because it would be crazy to try to come
  // up with names for these in a parameter list...
  // -- temperature dependence of density
  const double ka_, kb_, kc_, kd_;
  const double kT0_;

  // -- pressure dependence of density
  const double kalpha_, kp0_;

  static Utils::RegisteredFactoryWithState<EOS,EOSWater> factory_;

};

} // namespace
} // namespace
} // namespace

#endif
