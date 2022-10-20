/*
  Flow PK 

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

  Authors: Neil Carlson (nnc@lanl.gov), 
           Konstantin Lipnikov (lipnikov@lanl.gov)
           Daniil Svyatskiy (dasvyat@lanl.gov)

  Enforcing PDE properties on the discrete solution.
*/

#include <string>
#include <vector>

#include "Richards_PK.hh"

namespace Amanzi {
namespace Flow {

/* ******************************************************************
*
****************************************************************** */
void Richards_PK::CalculateCNLSLimiter_(
    const CompositeVector& wc, const CompositeVector& dwc_dp, double tol)
{ 
  auto& limiter = *cnls_limiter_->ViewComponent("cell");
  const auto& wcc = *wc.ViewComponent("cell");
  const auto& dwc_dpc = *dwc_dp.ViewComponent("cell");
  const auto& por = *S_->Get<CompositeVector>(porosity_key_).ViewComponent("cell");

  double alpha0, alpha1, wc_max, wc_min(0.0);
  double eps = tol * atm_pressure_; 
  for (int c = 0; c < ncells_owned; ++c) {
    wc_max = molar_rho_ * por[0][c];
    if (dwc_dpc[0][c] > 1e-16) {
      alpha0 = (wcc[0][c] - wc_min) / (dwc_dpc[0][c] * eps);
      alpha1 = (wc_max - wcc[0][c]) / (dwc_dpc[0][c] * eps); 
      limiter[0][c] = std::min(1.0, std::min(alpha0, alpha1));
    } else {
      limiter[0][c] = 1.0;
    }
  }
}


/* ******************************************************************
* 
****************************************************************** */
void Richards_PK::ApplyCNLSLimiter_()
{
}

}  // namespace Flow
}  // namespace Amanzi



