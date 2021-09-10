/*
  Shallow Water PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Rusanov model for numerical flux.
*/

#ifndef SHALLOW_WATER_NUMERICAL_FLUX_RUSANOV_HH_
#define SHALLOW_WATER_NUMERICAL_FLUX_RUSANOV_HH_

#include "Factory.hh"

#include "NumericalFlux.hh"

namespace Amanzi {
namespace ShallowWater {

class NumericalFlux_Rusanov : public NumericalFlux {
 public:
  explicit NumericalFlux_Rusanov(Teuchos::ParameterList& plist);
  ~NumericalFlux_Rusanov() {};
  
  virtual std::vector<double> Compute(
          const std::vector<double>& UL, const std::vector<double>& UR);

 private:
  static Utils::RegisteredFactory<NumericalFlux, NumericalFlux_Rusanov> factory_;
};

}  // namespace ShallowWater
}  // namespace Amanzi
 
#endif
