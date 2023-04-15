/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  Flow PK

  The Kozeny-Carman permeability factor depending on porosity.
*/

#ifndef AMANZI_FLOW_PERMEABILITY_KOZENY_CARMAN_HH_
#define AMANZI_FLOW_PERMEABILITY_KOZENY_CARMAN_HH_

#include <cmath>

#include "Teuchos_ParameterList.hpp"

#include "Permeability.hh"

namespace Amanzi {
namespace Flow {

class Permeability_KozenyCarman : public Permeability {
 public:
  explicit Permeability_KozenyCarman(Teuchos::ParameterList& plist)
  {
    phi_ref_ = plist.get<double>("undeformed soil porosity");
    factor_ = std::pow(1.0 - phi_ref_, 2) / std::pow(phi_ref_, 3);
  }
  ~Permeability_KozenyCarman(){};

  // required methods from the base class
  double Factor(double phi)
  {
    return factor_ * phi * phi * phi / (1.0 - phi) / (1.0 - phi);
  }

  double dFactordPorosity(double phi)
  {
    double a = phi * phi;
    double b = a * phi;
    return factor_ * (3 * a - b) / std::pow(1.0 - phi, 3);
  }

 private:
  double phi_ref_, factor_;
};

} // namespace Flow
} // namespace Amanzi

#endif
