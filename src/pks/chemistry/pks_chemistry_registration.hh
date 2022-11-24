/*
  Chemistry PK

  Copyright 2010-201x held jointly by LANS/LANL, LBNL, and PNNL. 
  Amanzi is released under the three-clause BSD License. 
  The terms of use and "as is" disclaimer for this license are 
  provided in the top-level COPYRIGHT file.

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
