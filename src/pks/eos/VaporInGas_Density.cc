/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Ethan Coon (coonet@ornl.gov)
*/

/*
  EOS

  EOS for a combination of air and vapor pressure. Mass density is not
  available, not because it can't be calculated, but because it depends
  upon omega. It's not really needed, and if it were, would not fit the
*/

#include "EOS_Density.hh"
#include "VaporInGas_Density.hh"
#include "EOSFactory.hh"

namespace Amanzi {
namespace AmanziEOS {

VaporInGas_Density::VaporInGas_Density(Teuchos::ParameterList& plist)
  : EOS_Density(plist.sublist("EOS parameters"))
{
  Teuchos::ParameterList gas_plist = plist.sublist("EOS parameters");
  EOSFactory<EOS_Density> eos_factory;
  gas_eos_ = eos_factory.Create(gas_plist);
}

} // namespace AmanziEOS
} // namespace Amanzi
