/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  EOS

  EOS for tabulated density for temperature between 0.5 and 800 C and
  pressure between 634 Pa and 110 MPa
*/

#include "DensityTabular.hh"
#include "LookupTable_Amanzi.hh"
#include "LookupTable_FEHM.hh"

namespace Amanzi {
namespace AmanziEOS {

DensityTabular::DensityTabular(Teuchos::ParameterList& plist) : EOS_Density(plist)
{
  std::string format = plist.get<std::string>("format", "Amanzi");
  if (format == "Amanzi")
    table_ = Teuchos::rcp(new LookupTable_Amanzi(plist));
  else 
    table_ = Teuchos::rcp(new LookupTable_FEHM(plist));
}

} // namespace AmanziEOS
} // namespace Amanzi
