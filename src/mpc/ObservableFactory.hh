/*
  Copyright 2010-202x held jointly by participating institutions.
  Amanzi is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Daniil Svyatskiy(dasvyat@lanl.gov)
*/

/*
  Multi-Process Coordinator

  A factory for creating observations.
*/

#ifndef AMANZI_OBSERVABLE_FACTORY_HH_
#define AMANZI_OBSERVABLE_FACTORY_HH_

#include "Epetra_MpiComm.h"
#include "Teuchos_ParameterList.hpp"

namespace Amanzi {

class Observable;

Teuchos::RCP<Observable>
CreateObservable(Teuchos::ParameterList& coord_plist,
                 Teuchos::ParameterList& observable_plist,
                 Teuchos::ParameterList& units_plist,
                 Teuchos::RCP<const AmanziMesh::Mesh> mesh);

} // namespace Amanzi

#endif
