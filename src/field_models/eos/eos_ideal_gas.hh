/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */

/*
  ATS

  EOS for an ideal gas

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

#ifndef FLOWRELATIONS_EOS_IDEAL_GAS_HH_
#define FLOWRELATIONS_EOS_IDEAL_GAS_HH_

#include "Teuchos_ParameterList.hpp"

#include "factory.hh"
#include "eos.hh"

namespace Amanzi {
namespace Flow {
namespace FlowRelations {

// Equation of State model
class EOSIdealGas : public EOS {

public:
  explicit
  EOSIdealGas(Teuchos::ParameterList& eos_plist);
  EOSIdealGas(const EOSIdealGas& other);

  virtual Teuchos::RCP<FieldModel> Clone() const;

  virtual double Density(double T, double p);
  virtual double DDensityDT(double T, double p);
  virtual double DDensityDp(double T, double p);

  virtual double molar_mass() const { return M_; }
  virtual bool is_molar_basis() const { return true; }

protected:
  virtual void InitializeFromPlist_();

  double R_;
  double M_;

private:
  static Utils::RegisteredFactory<FieldModel,EOSIdealGas> factory_;
};

} // namespace
} // namespace
} // namespace

#endif
