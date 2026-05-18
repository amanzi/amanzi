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

namespace Amanzi {
namespace AmanziEOS {

inline Teuchos::RCP<IAPWS95>
CreateIAPWS95(Teuchos::ParameterList& plist)
{
  if (plist.isParameter("csv table name")) {
    return Teuchos::rcp(new IAPWS95_Spline(plist));
  }
  return Teuchos::rcp(new IAPWS95(plist));
}

} // namespace AmanziEOS
} // namespace Amanzi

#endif
