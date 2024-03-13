/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Konstantin Lipnikov, Ethan Coon
*/

/*
  Process Kernels

  This is a base virtual class for process kernels. All physical
  kernels and MPCs must implement this interface for use within
  weak and strongly coupled hierarchies.
*/

#ifndef PK_IMPLICIT_HH_
#define PK_IMPLICIT_HH_

#include "Teuchos_RCP.hpp"

#include "ImplicitFn.hh"
#include "PK.hh"

namespace Amanzi {

class PK_Implicit : public PK, public ImplicitFn {};

} // namespace Amanzi

#endif
