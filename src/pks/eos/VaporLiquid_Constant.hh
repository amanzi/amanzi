/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  Constant vapor-liquid distribution coefficient
   
   kD = xg / xl 

  where xg and xl are gas and liquid more fractions in equilibrium.
*/

#ifndef AMANZI_EOS_VAPOR_LIQUID_CONSTANT_HH_
#define AMANZI_EOS_VAPOR_LIQUID_CONSTANT_HH_

#include <map>

#include "dbc.hh"
#include "Factory.hh"

#include "VaporLiquid.hh"

namespace Amanzi {
namespace AmanziEOS {

class VaporLiquid_Constant : public VaporLiquid {
 public:
  VaporLiquid_Constant(double kD) : kD_(kD){};

  virtual double k(double T) const { return kD_; }

  virtual double DkDT(double T) const { return 0.0; }

 private:
  double kD_;
};

} // namespace AmanziEOS
} // namespace Amanzi

#endif
