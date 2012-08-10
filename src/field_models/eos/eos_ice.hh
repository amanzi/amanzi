/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  ATS

  EOS for liquid ice.  See the permafrost physical properties notes for
  references and documentation of this EOS at:

  http://software.lanl.gov/ats/trac

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef FLOWRELATIONS_EOS_ICE_HH_
#define FLOWRELATIONS_EOS_ICE_HH_

#include "Teuchos_ParameterList.hpp"

#include "factory.hh"
#include "eos.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

// Equation of State model
class EOSIce : public EOS {

public:
  explicit
  EOSIce(Teuchos::ParameterList& eos_plist);

  EOSIce(const EOSIce& other);

  virtual Teuchos::RCP<FieldModel> Clone() const;

  virtual double Density(double T, double p);
  virtual double DDensityDT(double T, double p);
  virtual double DDensityDp(double T, double p);

  virtual double molar_mass() const { return M_; }
  virtual bool is_molar_basis() const { return true; }

private:
  virtual void InitializeFromPlist_();

  Teuchos::ParameterList eos_plist_;
  double M_;

  // constants for ice, hard-coded because it would be crazy to try to come
  // up with names for these in a parameter list...
  // -- temperature dependence of density
  const double ka_, kb_, kc_;
  const double kT0_;

  // -- pressure dependence of density
  const double kalpha_, kp0_;

  static Utils::RegisteredFactory<FieldModel,EOSIce> factory_;

};

} // namespace
} // namespace
} // namespace

#endif
