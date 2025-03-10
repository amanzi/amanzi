/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*!

Fugacity for the ideal gas is 1.

*/

#ifndef AMANZI_FUGACITY_IDEAL_GAS_HH_
#define AMANZI_FUGACITY_IDEAL_GAS_HH_

#include "Fugacity.hh"

namespace Amanzi {
namespace Multiphase {

class Fugacity_IdealGas : public Fugacity {
 public:
  Fugacity_IdealGas(const Teuchos::ParameterList& plist){};
  ~Fugacity_IdealGas() {};

  virtual double Value(double T) override { return 1.0; }
};

} // namespace Multiphase
} // namespace Amanzi

#endif
