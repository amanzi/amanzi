/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov (lipnikov@lanl.gov)
*/

/*
  EOS

  EOS for liquid water for temperature between 0.5 and 800 C and
  pressure between 634 Pa and 110 MPa
*/

#include "H2O_ViscosityTabular.hh"

namespace Amanzi {
namespace AmanziEOS {

H2O_ViscosityTabular::H2O_ViscosityTabular(Teuchos::ParameterList& eos_plist)
  : EOS_Viscosity(eos_plist)
{
  table_ = Teuchos::rcp(new LookupTable(eos_plist_));
}

} // namespace AmanziEOS
} // namespace Amanzi
