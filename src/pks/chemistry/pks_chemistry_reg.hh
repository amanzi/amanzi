/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors:
*/

/*
  Chemistry PK

  Self-registering factory for chemistry PKs.
*/

#ifdef ALQUIMIA_ENABLED
#  include "Alquimia_PK.hh"
#endif
#include "Amanzi_PK.hh"

namespace Amanzi {
namespace AmanziChemistry {

#ifdef ALQUIMIA_ENABLED
RegisteredPKFactory<Alquimia_PK> Alquimia_PK::reg_("chemistry alquimia");
#endif
RegisteredPKFactory<Amanzi_PK> Amanzi_PK::reg_("chemistry amanzi");

} // namespace AmanziChemistry
} // namespace Amanzi
