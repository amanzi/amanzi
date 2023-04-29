/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (ecoon@lanl.gov)
*/

/*
  Base class for vapor-liquid distribution coefficeint
   
   kD = xg / xl 

  where xg and xl are gas and liquid more fractions in equilibrium.
*/

#ifndef AMANZI_EOS_VAPOR_LIQUID_HH_
#define AMANZI_EOS_VAPOR_LIQUID_HH_

#include <map>

#include "dbc.hh"

#include "Factory.hh"

namespace Amanzi {
namespace AmanziEOS {

class VaporLiquid {
 public:
  VaporLiquid(){};
  virtual ~VaporLiquid(){};

  virtual double k(double T) const = 0;

  virtual double DkDT(double T) const = 0;
};

} // namespace AmanziEOS
} // namespace Amanzi

#endif
