// -----------------------------------------------------------------------------
// ATS
//
// License: see $ATS_DIR/COPYRIGHT
//
// Factory for taking coefficients for div-grad operators from cells to
// faces.
// -----------------------------------------------------------------------------

#include "errors.hh"

#include "UpwindFluxFactory.hh"
#include "upwind_total_flux.hh"
#include "upwind_flux_harmonic_mean.hh"
#include "upwind_flux_split_denominator.hh"
#include "upwind_flux_fo_cont.hh"
#include "upwind_cell_centered.hh"

namespace Amanzi {
namespace Operators {

Teuchos::RCP<Upwinding> UpwindFluxFactory::Create(Teuchos::ParameterList& oplist,
                                                  std::string pkname,
                                                  std::string cell_coef,
                                                  std::string face_coef,
                                                  std::string flux) {
  AMANZI_ASSERT(oplist.isSublist("overland conductivity model"));
  Teuchos::ParameterList sublist = oplist.sublist("overland conductivity model");
  std::string model_type = sublist.get<std::string>("overland conductivity type", "manning");
  double flux_eps = sublist.get<double>("upwind flux epsilon", 1.e-8);

  if (model_type == "manning") {
    return Teuchos::rcp(new UpwindTotalFlux(pkname, cell_coef, face_coef, flux, flux_eps));

  } else if (model_type == "manning harmonic mean") {
    return Teuchos::rcp(new UpwindFluxHarmonicMean(pkname, cell_coef, face_coef, flux, flux_eps));

  } else if (model_type == "manning split denominator") {
    std::string slope = oplist.get<std::string>("slope key", "slope_magnitude");
    std::string manning_coef = oplist.get<std::string>("coefficient key", "manning_coefficient");
    double slope_regularization = sublist.get<double>("slope regularization epsilon", 1.e-8);
    std::string ponded_depth = oplist.get<std::string>("height key", "ponded_depth");
    return Teuchos::rcp(new UpwindFluxSplitDenominator(pkname, cell_coef, face_coef, flux, flux_eps, slope, manning_coef, slope_regularization, ponded_depth));

  } else if (model_type == "manning ponded depth passthrough") {
    std::string slope = oplist.get<std::string>("slope key", "slope_magnitude");
    std::string manning_coef = oplist.get<std::string>("coefficient key", "manning_coefficient");
    double slope_regularization = sublist.get<double>("slope regularization epsilon", 1.e-8);
    double manning_exp = sublist.get<double>("Manning exponent");
    return Teuchos::rcp(new UpwindFluxFOCont(pkname, cell_coef, face_coef, flux, slope, manning_coef, "elevation", slope_regularization, manning_exp));

  } else if (model_type == "manning cell centered") {
    return Teuchos::rcp(new UpwindCellCentered(pkname, cell_coef, face_coef));
  } else {
    AMANZI_ASSERT(0);
    return Teuchos::null;
  }
}
  
}  // namespace Operators
}  // namespace Amanzi


