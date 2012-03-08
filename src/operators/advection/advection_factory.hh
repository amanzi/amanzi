/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
   ATS

   License: see $ATS_DIR/COPYRIGHT
   Author: Ethan Coon

   Interface for a general-purpose advection factory.
   ------------------------------------------------------------------------- */

#ifndef OPERATOR_ADVECTION_ADVECTION_FACTORY_HH_
#define OPERATOR_ADVECTION_ADVECTION_FACTORY_HH_

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Mesh.hh"
#include "composite_vector.hh"

#include "advection.hh"

namespace Amanzi {
namespace Operators {

class AdvectionFactory {

public:

  Teuchos::RCP<Advection> create(Teuchos::ParameterList& advect_plist,
          Teuchos::RCP<AmanziMesh::Mesh> mesh);
};

} // namespace Operators
} // namespace Amanzi

#endif
