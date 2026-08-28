/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  EOS

  Factory for IAPWS95 and derived models.
*/

#ifndef AMANZI_IAPWS95_FACTORY_HH_
#define AMANZI_IAPWS95_FACTORY_HH_

#include "Teuchos_ParameterList.hpp"

#include "IAPWS95.hh"
#include "IAPWS95_Spline.hh"
#include "IAPWS95_RaggedSplineRhoT.hh"

namespace Amanzi {
namespace AmanziEOS {

inline Teuchos::RCP<IAPWS95>
CreateIAPWS95(Teuchos::ParameterList& plist)
{
  if (plist.isParameter("csv table name")) {
    return Teuchos::rcp(new IAPWS95_Spline(plist));
  } else if (plist.isParameter("use ragged spline")) {
    IAPWS95_RaggedSplineRhoT::Options opt;
    opt.T_min = 280.0;
    opt.T_max = 950.0;
    opt.rho_min = 0.5;
    opt.rho_max = 1100.0;
    opt.initial_rho_intervals = 60;
    opt.initial_T_intervals = 70;
    opt.max_rho_intervals = 130;
    opt.max_T_intervals = 180;
    opt.extension_cells = 4.0;
    opt.anisotropic_refinement_fraction = 0.25;

    Teuchos::ParameterList plist;
    auto spline = Teuchos::rcp(new IAPWS95_RaggedSplineRhoT(plist, opt));
    spline->InitializeSharedData();

    return spline;
  } else {
    return Teuchos::rcp(new IAPWS95(plist));
  }
}

} // namespace AmanziEOS
} // namespace Amanzi

#endif
