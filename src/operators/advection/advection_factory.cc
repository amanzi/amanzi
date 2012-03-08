/* -*-  mode: c++; c-default-style: "google"; indent-tabs-mode: nil -*- */
/* -------------------------------------------------------------------------
   ATS

   License: see $ATS_DIR/COPYRIGHT
   Author: Ethan Coon

   Interface for a general-purpose advection factory.
   ------------------------------------------------------------------------- */

#include "advection_donor_upwind.hh"
#include "advection_factory.hh"

namespace Amanzi {
namespace Operators {

Teuchos::RCP<Advection> AdvectionFactory::create(Teuchos::ParameterList& plist,
        Teuchos::RCP<AmanziMesh::Mesh> mesh) {
  std::string method = plist.get<string>("Advection method");

  if (method == "donor upwind") {
    return Teuchos::rcp(new AdvectionDonorUpwind(plist, mesh));
  }
};

} // namespace Operators
} // namespace Amanzi
