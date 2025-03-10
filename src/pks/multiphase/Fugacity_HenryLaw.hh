/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*!

Fugacity model for amount of dissolved gas in a liquid phase.

*/

#ifndef AMANZI_FUGACITY_HENRY_LAW_HH_
#define AMANZI_FUGACITY_HENRY_LAW_HH_

#include "Fugacity.hh"

namespace Amanzi {
namespace Multiphase {

class Fugacity_HenryLaw : public Fugacity {
 public:
  Fugacity_HenryLaw(const Teuchos::ParameterList& plist) {
    kH_ = plist.get<double>("Henry constant");
  }
  ~Fugacity_HenryLaw() {};

  virtual double Value(double T) override { return kH_; }

 private:
  double kH_;
};

} // namespace Multiphase
} // namespace Amanzi

#endif
