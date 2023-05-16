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

namespace Amanzi {
namespace AmanziEOS {

DensityTabular::DensityTabular(Teuchos::ParameterList& eos_plist) : EOS_Density(eos_plist)
{
  table_ = Teuchos::rcp(new LookupTable(eos_plist_));
}

} // namespace AmanziEOS
} // namespace Amanzi
