/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon
*/

/*
  Energy

  EOS for liquid water based on the International Association for the 
  Properties of Water and Steam (IAPWS), Industrial Formulation 1997.
*/

#ifndef AMANZI_ENERGY_IEM_IAPWS97_HH_
#define AMANZI_ENERGY_IEM_IAPWS97_HH_

#include "Teuchos_ParameterList.hpp"

#include "IAPWS97.hh"

#include "IEM.hh"
#include "Factory.hh"

namespace Amanzi {
namespace Energy {

class IEM_IAPWS97 : public IEM {
 public:
  explicit IEM_IAPWS97(Teuchos::ParameterList& plist);

  virtual double InternalEnergy(double T, double p) override;
  virtual double DInternalEnergyDT(double T, double p) override;
  virtual double DInternalEnergyDp(double T, double p) override;

 private:
  Teuchos::RCP<AmanziEOS::IAPWS97> eos_;
  static Utils::RegisteredFactory<IEM, IEM_IAPWS97> reg_;
};

} // namespace Energy
} // namespace Amanzi

#endif
