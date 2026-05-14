/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  EOS

  Revised Release on the IAPWS-95 formulation
  for the Thermodynamic Properties of Water and Steam.
*/

#ifndef AMANZI_IAPWS95_SPLINE_HH_
#define AMANZI_IAPWS95_SPLINE_HH_

#include <tuple>

#include "Teuchos_ParameterList.hpp"

// Amanzi
#include "SplineCubicNotAKnot2D.hh"

// Amanzi::EOS
#include "IAPWS95.hh"

namespace Amanzi {
namespace AmanziEOS {

class IAPWS95_Spline : public IAPWS95 {
 public:
  IAPWS95_Spline(Teuchos::ParameterList& plist);
  ~IAPWS95_Spline() {};

  virtual std::array<double, 6> ResidualPart(double rho, double T) override;

 private:
  double Build_(double tau, double delta) const;

  void ReadTable_(const std::string& filename);
  std::vector<std::string> SplitLine_(const std::string& line);
  bool CheckStrictMonotonicity_(const std::vector<double>& x);

 private:
  std::vector<double> tau_;
  std::vector<double> delta_;
  std::vector<std::vector<double>> values_;

  WhetStone::SplineCubicNotAKnot2D spline_;
  int nx_, ny_;
};

} // namespace AmanziEOS
} // namespace Amanzi

#endif


